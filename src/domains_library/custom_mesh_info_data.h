#ifndef __CUSTOM_MESH_INFO_DATA_H
#define __CUSTOM_MESH_INFO_DATA_H

#include <stdbool.h>
#include "../common_types/common_types.h"

//Written by Leto Riebel

struct paci_patch_mesh_info {
	int is_paci;
	int layer;
  int infarct_zone;
  float apicobasal;
	int fast_endo;
};

#define PACI_INFO(grid_cell) (struct paci_patch_mesh_info *)grid_cell->mesh_extra_info
//IS_PACI: 1 = ventricular-like, 2 = atrial-like, 3 = nodal-like
#define IS_PACI(grid_cell) (PACI_INFO(grid_cell))->is_paci
//LAYER: 1 = endocardial, 2 = midmyocardial (not use dhere), 3 = epicardial
#define LAYER(grid_cell) (PACI_INFO(grid_cell))->layer
//INFARCT_ZONE: 0 = healthy, 1 = infarct center, 2 = border zone
#define INFARCT_ZONE(grid_cell) (PACI_INFO(grid_cell))->infarct_zone
//APICOBASAL: float between 0 and 1 which in the model file is normalised to produce IKs scaling factor between 0.2 and 5.0
#define APICOBASAL(grid_cell) (PACI_INFO(grid_cell))->apicobasal
//FAST_ENDO: 0 = myocardium, 1= fast conducting endocardial layer (assigned but not used here)
#define FAST_ENDO(grid_cell) (PACI_INFO(grid_cell))->fast_endo

#define INITIALIZE_PACI_INFO(grid_cell)                                                                            \
    do {                                                                                                           \
        size_t __size__ = sizeof (struct paci_patch_mesh_info);                                                    \
        (grid_cell)->mesh_extra_info = malloc (__size__);                                                          \
        (grid_cell)->mesh_extra_info_size = __size__;                                                              \
        IS_PACI ((grid_cell)) = 4;                                                                                 \
        LAYER ((grid_cell)) = 4;                                                                             \
        INFARCT_ZONE ((grid_cell)) = 3;                                                                                   \
        APICOBASAL ((grid_cell)) = 1;                                                                                   \
        FAST_ENDO ((grid_cell)) = 2;                                                                               \
} while (0)

#endif /* __MESH_INFO_DATA_H */
