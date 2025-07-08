#ifndef CUSTOM_SCHEDULER_H
#define CUSTOM_SCHEDULER_H

#include <stdbool.h>
#include "../oai/openair2/LAYER2/NR_MAC_gNB/nr_mac_gNB.h"

typedef struct UEsched_s {
  float coef;
  NR_UE_info_t * UE;
} UEsched_t;

extern UEsched_t UEsched_list[64];

extern bool use_custom_scheduler;

#endif