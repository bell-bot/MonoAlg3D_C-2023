#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../alg/grid/grid.h"
#include "../config/assembly_matrix_config.h"
#include "../libraries_common/common_data_structures.h"
#include "../utils/file_utils.h"
#include "../utils/utils.h"
#include "../domains_library/custom_mesh_info_data.h"
#include "purkinje_coupling_matrix_assembly.c"

#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#endif

//Adapted by Leto Riebel

//Use if no Purkinje system
ASSEMBLY_MATRIX(stem_cell_coupling_assembly_matrix) {
    // [TISSUE]
    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    initialize_diagonal_elements(the_solver, the_grid);

    //      D tensor    //
    // | sx    sxy   sxz |
    // | sxy   sy    syz |
    // | sxz   syz   sz  |
    real_cpu D[3][3];
    int i;

    real_cpu healthy_sigma_l = 0.0;
    real_cpu healthy_sigma_t = 0.0;
    real_cpu healthy_sigma_n = 0.0;
    real_cpu stem_cell_sigma_l = 0.0;
    real_cpu stem_cell_sigma_t = 0.0;
    real_cpu stem_cell_sigma_n = 0.0;
    real_cpu infarction_sigma_l = 0.0;
    real_cpu infarction_sigma_t = 0.0;
    real_cpu infarction_sigma_n = 0.0;
    real_cpu fast_endo_scale = 0.0;

    char *fiber_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(fiber_file, config, "fibers_file");

    bool fibers_in_mesh = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(fibers_in_mesh, config, "fibers_in_mesh");
    struct fiber_coords *fibers = NULL;

    uint32_t infarct_stage = 4;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(uint32_t, infarct_stage, config, "infarct_stage");
    //Healthy ventricular conductivities
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, healthy_sigma_l, config, "sigma_l");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, healthy_sigma_t, config, "sigma_t");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, healthy_sigma_n, config, "sigma_n");
    //Stem cell conductivities
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, stem_cell_sigma_l, config, "stem_cell_sigma_l");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, stem_cell_sigma_t, config, "stem_cell_sigma_t");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, stem_cell_sigma_n, config, "stem_cell_sigma_n");
    //Infarction conductivities
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, infarction_sigma_l, config, "infarction_sigma_l");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, infarction_sigma_t, config, "infarction_sigma_t");
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, infarction_sigma_n, config, "infarction_sigma_n");
    //Fast endocardial layer - scale factor, set to 1 if using Purkinje
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, fast_endo_scale, config, "fast_endo_scale");

    real_cpu *f = NULL;
    real_cpu *s = NULL;
    real_cpu *n = NULL;

    if(fiber_file) {
        log_info("Loading mesh fibers\n");
        fibers = read_fibers(fiber_file, false);
    } else if(!fibers_in_mesh) {
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(f, config, "f", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(s, config, "s", 3);
        GET_PARAMETER_VECTOR_VALUE_OR_USE_DEFAULT(n, config, "n", 3);

        if(!f) {
            f = malloc(sizeof(real_cpu)*3);
            f[0] = 1.0;
            f[1] = 0.0;
            f[2] = 0.0;
        }

        if(!s) {
            s = malloc(sizeof(real_cpu)*3);
            s[0] = 0.0;
            s[1] = 1.0;
            s[2] = 0.0;
        }

        if(!n) {
            n = malloc(sizeof(real_cpu)*3);
            n[0] = 0.0;
            n[1] = 0.0;
            n[2] = 1.0;
        }
    }

    OMP(parallel for private(D))
    for(i = 0; i < num_active_cells; i++) {
        real_cpu sigma_l = 0.0;
        real_cpu sigma_t = 0.0;
        real_cpu sigma_n = 0.0;
        if((IS_PACI(ac[i]) > 0) && (IS_PACI(ac[i]) < 4)) {
            sigma_l = stem_cell_sigma_l;
            sigma_t = stem_cell_sigma_t;
            sigma_n = stem_cell_sigma_n;
        } else if((INFARCT_ZONE(ac[i]) >= 3) || (IS_PACI(ac[i]) >= 4) || (FAST_ENDO(ac[i]) >= 2) || (infarct_stage >= 4) || (LAYER(ac[i]) >= 4)) {
            log_error_and_exit("Infarction zone or is_paci has not been set for cell in index %d - %lf, %lf, %lf\n", i, ac[i]->center.x, ac[i]->center.y, ac[i]->center.z);
        } else if((infarct_stage < 1) || (INFARCT_ZONE(ac[i]) == 0)) {
            sigma_l = healthy_sigma_l;
            sigma_t = healthy_sigma_t;
            sigma_n = healthy_sigma_n;
        } else if((INFARCT_ZONE(ac[i]) > 0) && (infarct_stage > 0)) {
            sigma_l = infarction_sigma_l;
            sigma_t = infarction_sigma_t;
            sigma_n = infarction_sigma_n;
        } else {
            log_error_and_exit("Error in assigning conductivities, mismatch in infarct zone or infarct stage.");
        }

        if(fibers) {
            int fiber_index = ac[i]->original_position_in_file;

            if(fiber_index == -1) {
                log_error_and_exit("fiber_index should not be -1, but it is for cell in index %d - %lf, %lf, %lf\n", i, ac[i]->center.x, ac[i]->center.y, ac[i]->center.z);
            }

            if(sigma_t == sigma_n) {
                calc_tensor2(D, fibers[fiber_index].f, sigma_l, sigma_t);
            }
            else {
                calc_tensor(D, fibers[fiber_index].f, fibers[fiber_index].s, fibers[fiber_index].n, sigma_l, sigma_t, sigma_n);
            }
            ac[i]->sigma.fibers = fibers[fiber_index];
        }
        else if(fibers_in_mesh) {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, ac[i]->sigma.fibers.f, sigma_l, sigma_t);
            }
            else {
                calc_tensor(D, ac[i]->sigma.fibers.f, ac[i]->sigma.fibers.s, ac[i]->sigma.fibers.n, sigma_l, sigma_t, sigma_n);
            }

        }
        else {
            if(sigma_t == sigma_n) {
                calc_tensor2(D, f, sigma_l, sigma_t);
            }
            else {
                calc_tensor(D, f, s, n, sigma_l, sigma_t, sigma_n);
            }
        }

        if ((FAST_ENDO(ac[i]) == 1) && (IS_PACI(ac[i]) == 0)) {
            ac[i]->sigma.x = D[0][0] * fast_endo_scale;
            ac[i]->sigma.y = D[1][1] * fast_endo_scale;
            ac[i]->sigma.z = D[2][2] * fast_endo_scale;

            ac[i]->sigma.xy = D[0][1] * fast_endo_scale;
            ac[i]->sigma.xz = D[0][2] * fast_endo_scale;
            ac[i]->sigma.yz = D[1][2] * fast_endo_scale;
        }
        // Normal type of cell
        else {
            ac[i]->sigma.x = D[0][0];
            ac[i]->sigma.y = D[1][1];
            ac[i]->sigma.z = D[2][2];

            ac[i]->sigma.xy = D[0][1];
            ac[i]->sigma.xz = D[0][2];
            ac[i]->sigma.yz = D[1][2];
        }
    }

    OMP(parallel for)
    for(i = 0; i < num_active_cells; i++) {
        fill_discretization_matrix_elements_aniso(ac[i]);
    }

    free(f);
    free(s);
    free(n);
}
