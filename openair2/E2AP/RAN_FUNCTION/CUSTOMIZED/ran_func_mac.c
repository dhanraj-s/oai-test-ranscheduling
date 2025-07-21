/*
 * Licensed to the OpenAirInterface (OAI) Software Alliance under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The OpenAirInterface Software Alliance licenses this file to You under
 * the OAI Public License, Version 1.1  (the "License"); you may not use this file
 * except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.openairinterface.org/?page_id=698
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *-------------------------------------------------------------------------------
 * For more information about the OpenAirInterface (OAI) Software Alliance:
 *      contact@openairinterface.org
 */ 
#include <unistd.h>
#include "ran_func_mac.h"
#include <assert.h>
#include <stdio.h>
#include <../../../../openair1/PHY/defs_RU.h>
#include <../../../../executables/softmodem-common.h>
#include <../../../LAYER2/nr_rlc/nr_rlc_oai_api.h>
#include "../../../NR_PHY_INTERFACE/NR_IF_Module.h"
#include "../../../../executables/custom_scheduler.h"

static
const int mod_id = 0;


bool read_mac_sm(void* data)
{
  //printf("\n\n\n\n\nIN THE REPORT CALLBACK!\n\n\n\n\n");
  assert(data != NULL);

  mac_ind_data_t* mac = (mac_ind_data_t*)data;
  //fill_mac_ind_data(mac);

  mac->msg.tstamp = time_now_us();

  NR_UEs_t *UE_info = &RC.nrmac[mod_id]->UE_info;

  //printf("timestamp_rx: %ld\n", RC.ru[0]->proc.timestamp_rx);
  //printf("frame_rx: %ld\n", RC.ru[0]->proc.frame_rx);
  //printf("tti_rx: %ld\n\n", RC.ru[0]->proc.tti_rx);

  size_t num_ues = 0;
  UE_iterator(UE_info->connected_ue_list, ue) {
    if (ue)
      num_ues += 1;
  }

  mac->msg.len_ue_stats = num_ues;
  if(mac->msg.len_ue_stats > 0){
    mac->msg.ue_stats = calloc(mac->msg.len_ue_stats, sizeof(mac_ue_stats_impl_t));
    assert(mac->msg.ue_stats != NULL && "Memory exhausted" );
  }

  size_t i = 0; //TODO
  UE_iterator(UE_info->connected_ue_list, UE) {
    const NR_UE_sched_ctrl_t* sched_ctrl = &UE->UE_sched_ctrl;
    mac_ue_stats_impl_t* rd = &mac->msg.ue_stats[i];

    rd->frame = RC.nrmac[mod_id]->frame;
    rd->slot = 0; // previously had slot info, but the gNB runs multiple slots
                  // in parallel, so this has no real meaning
    rd->dl_aggr_tbs = UE->mac_stats.dl.total_bytes;
    rd->ul_aggr_tbs = UE->mac_stats.ul.total_bytes;

    if (is_dl_slot(rd->slot, &RC.nrmac[mod_id]->frame_structure)) {
      rd->dl_curr_tbs = UE->mac_stats.dl.current_bytes;
      rd->dl_sched_rb = UE->mac_stats.dl.current_rbs;
    }
    if (is_ul_slot(rd->slot, &RC.nrmac[mod_id]->frame_structure)) {
      rd->ul_curr_tbs = UE->mac_stats.ul.current_bytes;
      rd->ul_sched_rb = sched_ctrl->sched_pusch.rbSize;
    }

    rd->rnti = UE->rnti;
    rd->dl_aggr_prb = UE->mac_stats.dl.total_rbs;
    rd->ul_aggr_prb = UE->mac_stats.ul.total_rbs;
    rd->dl_aggr_retx_prb = UE->mac_stats.dl.total_rbs_retx;
    rd->ul_aggr_retx_prb = UE->mac_stats.ul.total_rbs_retx;

    rd->dl_aggr_bytes_sdus = UE->mac_stats.dl.lc_bytes[3];
    rd->ul_aggr_bytes_sdus = UE->mac_stats.ul.lc_bytes[3];

    rd->dl_aggr_sdus = UE->mac_stats.dl.num_mac_sdu;
    rd->ul_aggr_sdus = UE->mac_stats.ul.num_mac_sdu;

    rd->pusch_snr = (float) sched_ctrl->pusch_snrx10 / 10; //: float = -64;
    rd->pucch_snr = (float) sched_ctrl->pucch_snrx10 / 10; //: float = -64;

    rd->wb_cqi = sched_ctrl->CSI_report.cri_ri_li_pmi_cqi_report.wb_cqi_1tb;
    rd->dl_mcs1 = sched_ctrl->dl_bler_stats.mcs;
    rd->dl_bler = sched_ctrl->dl_bler_stats.bler;
    rd->ul_mcs1 = sched_ctrl->ul_bler_stats.mcs;
    rd->ul_bler = sched_ctrl->ul_bler_stats.bler;
    rd->dl_mcs2 = 0;
    rd->ul_mcs2 = 0;
    rd->phr = sched_ctrl->ph;

    //printf("sched_ctrl->estimated_ul_buffer: %u\n", sched_ctrl->estimated_ul_buffer);
    //printf("sched_ctrl->sched_ul_bytes: %u\n", sched_ctrl->sched_ul_bytes);
    // Bug: difference of unsigned quantities when negative will give a very large unsigned number.
    const uint32_t bufferSize = sched_ctrl->estimated_ul_buffer - sched_ctrl->sched_ul_bytes;
    //printf("bufferSize: %u\n", bufferSize);
    rd->bsr = bufferSize;

    const size_t numDLHarq = 4;
    rd->dl_num_harq = numDLHarq;
    for (uint8_t j = 0; j < numDLHarq; ++j)
      rd->dl_harq[j] = UE->mac_stats.dl.rounds[j];
    rd->dl_harq[numDLHarq] = UE->mac_stats.dl.errors;

    const size_t numUlHarq = 4;
    rd->ul_num_harq = numUlHarq;
    for (uint8_t j = 0; j < numUlHarq; ++j)
      rd->ul_harq[j] = UE->mac_stats.ul.rounds[j];
    rd->ul_harq[numUlHarq] = UE->mac_stats.ul.errors;

    ++i;
  }

  return num_ues > 0;
}

void read_mac_setup_sm(void* data)
{
  assert(data != NULL);
  assert(0 !=0 && "Not supported");
}

void create_ue_sched_list( UEsched_t *UE_sched, mac_ctrl_msg_t mac_ctrl_msg, NR_UE_info_t **UE_list, pthread_mutex_t *mutex) {
  pthread_mutex_lock(mutex);
  for(int i=0; i<64; ++i) {
    UE_sched[i].UE = NULL;
    UE_sched[i].coef = 1.0;
  }

  for(int i=0; i<mac_ctrl_msg.num_users; ++i) {
    int cur_rnti = mac_ctrl_msg.resource_alloc[i].user_id;
    UE_iterator(UE_list, UE) {
      if(UE->rnti == cur_rnti) {
        UE_sched[i].UE = UE;
        UE_sched[i].coef = 1.0;
        break;
      }
    }
  }

  use_custom_scheduler = true;
  pthread_mutex_unlock(mutex);
}

int get_mcs_from_ctrl(int rnti, mac_ctrl_msg_t mac_ctrl_msg) {
  for(int i=0; i<mac_ctrl_msg.num_users; ++i) {
    if(mac_ctrl_msg.resource_alloc[i].user_id == rnti)
      return mac_ctrl_msg.resource_alloc[i].mcs;
  }
}

int get_rbSize_from_ctrl(int rnti, mac_ctrl_msg_t mac_ctrl_msg) {
  for(int i=0; i<mac_ctrl_msg.num_users; ++i) {
    if(mac_ctrl_msg.resource_alloc[i].user_id == rnti)
      return mac_ctrl_msg.resource_alloc[i].num_rb;
  }
}


sm_ag_if_ans_t write_ctrl_mac_sm(void const* data)
{
  assert(data != NULL);
  
  /*Get the preprocessor parameters ready*/
  frame_t frame = RC.ru[0]->proc.frame_tx; slot_t slot = RC.ru[0]->proc.tti_tx;
  module_id_t module_id = 0;

  /*Get the MAC Control message for the scheduling decisions*/
  mac_ctrl_req_data_t *mac_data = data;
  mac_ctrl_msg_t mac_ctrl_msg = mac_data->msg;

  /*Create the UE order for scheduling.*/
  gNB_MAC_INST *mac = RC.nrmac[mod_id];
  NR_UEs_t *UE_info = &mac->UE_info;

  pthread_mutex_t list_mutex;
  pthread_mutex_init(&list_mutex, NULL);
  create_ue_sched_list(UEsched_list, mac_ctrl_msg, UE_info->connected_ue_list, &list_mutex);
  pthread_mutex_destroy(&list_mutex);
  sm_ag_if_ans_t ans = {0};
  return ans;
}

