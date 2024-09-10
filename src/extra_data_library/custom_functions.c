#include <unistd.h>

#include "../config/extra_data_config.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../utils/file_utils.h"
#include "../domains_library/custom_mesh_info_data.h"

//Written by Leto Riebel

//Sets extra data array for stem cell region, transmurality, apicobasal gradient,
//and infarct zone from the mesh file and infarct stage from the .ini file
SET_EXTRA_DATA(set_extra_data_with_stem_cells) {
    uint32_t num_active_cells = the_grid->num_active_cells;
    struct cell_node ** ac = the_grid->active_cells;

    *extra_data_size = sizeof(real)*((num_active_cells)*5 + 45);
    real *extra_data = (real*)malloc(*extra_data_size);

    uint32_t infarct_stage = 5;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(uint32_t, infarct_stage, config, "infarct_stage");

    //These are applied to ALL cells including iPSC-CMs, useful for simulating drugs
    real INa_mult = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INa_mult, config, "INa_mult");
    real IKr_mult = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKr_mult, config, "IKr_mult");
    real ICaL_mult = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICaL_mult, config, "ICaL_mult");
    real INaL_mult = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaL_mult, config, "INaL_mult");
    real IKs_mult = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKs_mult, config, "IKs_mult");
    real Ito_mult = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ito_mult, config, "Ito_mult");
    real IK1_mult = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IK1_mult, config, "IK1_mult");

    //Theses are applied to the REMOTE ZONE, i.e. non-infarcted tissue and not including iPSC-CMs
    real INa_mult_RZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INa_mult_RZ, config, "INa_mult_RZ");
    real IKr_mult_RZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKr_mult_RZ, config, "IKr_mult_RZ");
    real ICaL_mult_RZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICaL_mult_RZ, config, "ICaL_mult_RZ");
    real INaL_mult_RZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaL_mult_RZ, config, "INaL_mult_RZ");
    real IKs_mult_RZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKs_mult_RZ, config, "IKs_mult_RZ");
    real Ito_mult_RZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ito_mult_RZ, config, "Ito_mult_RZ");
    real IK1_mult_RZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IK1_mult_RZ, config, "IK1_mult_RZ");
    real aCaMK_mult_RZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, aCaMK_mult_RZ, config, "aCaMK_mult_RZ");
    real tau_relp_mult_RZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, tau_relp_mult_RZ, config, "tau_relp_mult_RZ");
    real ICab_mult_RZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICab_mult_RZ, config, "ICab_mult_RZ");
    real Iup_mult_RZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Iup_mult_RZ, config, "Iup_mult_RZ");
    real IKCa_mult_RZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKCa_mult_RZ, config, "IKCa_mult_RZ");
    real IClCa_mult_RZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IClCa_mult_RZ, config, "IClCa_mult_RZ");

    //Theses are applied to the INFARCT CENTER, not including iPSC-CMs
    real INa_mult_IZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INa_mult_IZ, config, "INa_mult_IZ");
    real IKr_mult_IZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKr_mult_IZ, config, "IKr_mult_IZ");
    real ICaL_mult_IZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICaL_mult_IZ, config, "ICaL_mult_IZ");
    real INaL_mult_IZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaL_mult_IZ, config, "INaL_mult_IZ");
    real IKs_mult_IZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKs_mult_IZ, config, "IKs_mult_IZ");
    real Ito_mult_IZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ito_mult_IZ, config, "Ito_mult_IZ");
    real IK1_mult_IZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IK1_mult_IZ, config, "IK1_mult_IZ");
    real aCaMK_mult_IZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, aCaMK_mult_IZ, config, "aCaMK_mult_IZ");
    real tau_relp_mult_IZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, tau_relp_mult_IZ, config, "tau_relp_mult_IZ");
    real ICab_mult_IZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICab_mult_IZ, config, "ICab_mult_IZ");
    real Ko_mult_IZ = 5;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ko_mult_IZ, config, "Ko_mult_IZ");

    //Theses are applied to the BORDER ZONE, not including iPSC-CMs
    real INa_mult_BZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INa_mult_BZ, config, "INa_mult_BZ");
    real IKr_mult_BZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKr_mult_BZ, config, "IKr_mult_BZ");
    real ICaL_mult_BZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICaL_mult_BZ, config, "ICaL_mult_BZ");
    real INaL_mult_BZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, INaL_mult_BZ, config, "INaL_mult_BZ");
    real IKs_mult_BZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKs_mult_BZ, config, "IKs_mult_BZ");
    real Ito_mult_BZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Ito_mult_BZ, config, "Ito_mult_BZ");
    real IK1_mult_BZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IK1_mult_BZ, config, "IK1_mult_BZ");
    real aCaMK_mult_BZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, aCaMK_mult_BZ, config, "aCaMK_mult_BZ");
    real tau_relp_mult_BZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, tau_relp_mult_BZ, config, "tau_relp_mult_BZ");
    real ICab_mult_BZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, ICab_mult_BZ, config, "ICab_mult_BZ");
    real Iup_mult_BZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, Iup_mult_BZ, config, "Iup_mult_BZ");
    real IKCa_mult_BZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IKCa_mult_BZ, config, "IKCa_mult_BZ");
    real IClCa_mult_BZ = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real, IClCa_mult_BZ, config, "IClCa_mult_BZ");

    //Layer, infarct zone, apicobasal and fast-endo are set for each cell
    uint32_t i;
    OMP(parallel for)
    for (i = 0; i < num_active_cells; i++) {
        extra_data[i] = IS_PACI(ac[i]); //tor 0, paci 1, atrial 2, nodal 3
        extra_data[i+num_active_cells] = LAYER(ac[i]); //fast-endo 0, endo 1, mid 2, epi 3
        extra_data[i+(2*num_active_cells)] = INFARCT_ZONE(ac[i]); //healthy 0, infarct 1, border 2, remote 3
        extra_data[i+(3*num_active_cells)] = APICOBASAL(ac[i]); //float number between 0 and 1
        extra_data[i+(4*num_active_cells)] = FAST_ENDO(ac[i]); //float number between 0 and 1
    }
    //Remodellings are set globally and applied only to the relevant cells, identified by layer and infarct zone
    extra_data[(5*num_active_cells)] = infarct_stage;
    extra_data[(5*num_active_cells + 1)] = INa_mult;
    extra_data[(5*num_active_cells + 2)] = IKr_mult;
    extra_data[(5*num_active_cells + 3)] = ICaL_mult;
    extra_data[(5*num_active_cells + 4)] = INaL_mult;
    extra_data[(5*num_active_cells + 5)] = IKs_mult;
    extra_data[(5*num_active_cells + 6)] = Ito_mult;
    extra_data[(5*num_active_cells + 7)] = IK1_mult;

    extra_data[(5*num_active_cells + 8)] = INa_mult_RZ;
    extra_data[(5*num_active_cells + 9)] = IKr_mult_RZ;
    extra_data[(5*num_active_cells + 10)] = ICaL_mult_RZ;
    extra_data[(5*num_active_cells + 11)] = INaL_mult_RZ;
    extra_data[(5*num_active_cells + 12)] = IKs_mult_RZ;
    extra_data[(5*num_active_cells + 13)] = Ito_mult_RZ;
    extra_data[(5*num_active_cells + 14)] = IK1_mult_RZ;
    extra_data[(5*num_active_cells + 15)] = aCaMK_mult_RZ;
    extra_data[(5*num_active_cells + 16)] = tau_relp_mult_RZ;
    extra_data[(5*num_active_cells + 17)] = ICab_mult_RZ;
    extra_data[(5*num_active_cells + 18)] = Iup_mult_RZ;
    extra_data[(5*num_active_cells + 19)] = IKCa_mult_RZ;
    extra_data[(5*num_active_cells + 20)] = IClCa_mult_RZ;

    extra_data[(5*num_active_cells + 21)] = INa_mult_IZ;
    extra_data[(5*num_active_cells + 22)] = IKr_mult_IZ;
    extra_data[(5*num_active_cells + 23)] = ICaL_mult_IZ;
    extra_data[(5*num_active_cells + 24)] = INaL_mult_IZ;
    extra_data[(5*num_active_cells + 25)] = IKs_mult_IZ;
    extra_data[(5*num_active_cells + 26)] = Ito_mult_IZ;
    extra_data[(5*num_active_cells + 27)] = IK1_mult_IZ;
    extra_data[(5*num_active_cells + 28)] = aCaMK_mult_IZ;
    extra_data[(5*num_active_cells + 29)] = tau_relp_mult_IZ;
    extra_data[(5*num_active_cells + 30)] = ICab_mult_IZ;
    extra_data[(5*num_active_cells + 31)] = Ko_mult_IZ;

    extra_data[(5*num_active_cells + 32)] = INa_mult_BZ;
    extra_data[(5*num_active_cells + 33)] = IKr_mult_BZ;
    extra_data[(5*num_active_cells + 34)] = ICaL_mult_BZ;
    extra_data[(5*num_active_cells + 35)] = INaL_mult_BZ;
    extra_data[(5*num_active_cells + 36)] = IKs_mult_BZ;
    extra_data[(5*num_active_cells + 37)] = Ito_mult_BZ;
    extra_data[(5*num_active_cells + 38)] = IK1_mult_BZ;
    extra_data[(5*num_active_cells + 39)] = aCaMK_mult_BZ;
    extra_data[(5*num_active_cells + 40)] = tau_relp_mult_BZ;
    extra_data[(5*num_active_cells + 41)] = ICab_mult_BZ;
    extra_data[(5*num_active_cells + 42)] = Iup_mult_BZ;
    extra_data[(5*num_active_cells + 43)] = IKCa_mult_BZ;
    extra_data[(5*num_active_cells + 44)] = IClCa_mult_BZ;

    return extra_data;
}
