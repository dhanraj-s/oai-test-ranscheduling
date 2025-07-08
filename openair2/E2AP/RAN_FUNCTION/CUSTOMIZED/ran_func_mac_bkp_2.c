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

static void nr_store_dlsch_buffer(module_id_t module_id, frame_t frame, slot_t slot)
{
  UE_iterator(RC.nrmac[module_id]->UE_info.connected_ue_list, UE) {
    NR_UE_sched_ctrl_t *sched_ctrl = &UE->UE_sched_ctrl;
    sched_ctrl->num_total_bytes = 0;
    sched_ctrl->dl_pdus_total = 0;

    /* loop over all activated logical channels */
    // Note: DL_SCH_LCID_DCCH, DL_SCH_LCID_DCCH1, DL_SCH_LCID_DTCH
    for (int i = 0; i < seq_arr_size(&sched_ctrl->lc_config); ++i) {
      const nr_lc_config_t *c = seq_arr_at(&sched_ctrl->lc_config, i);
      const int lcid = c->lcid;
      const uint16_t rnti = UE->rnti;
      LOG_D(NR_MAC, "In %s: UE %x: LCID %d\n", __FUNCTION__, rnti, lcid);
      if (lcid == DL_SCH_LCID_DTCH && nr_timer_is_active(&sched_ctrl->transm_interrupt))
        continue;
      start_meas(&RC.nrmac[module_id]->rlc_status_ind);
      sched_ctrl->rlc_status[lcid] = nr_mac_rlc_status_ind(rnti, frame, lcid);
      stop_meas(&RC.nrmac[module_id]->rlc_status_ind);

      if (sched_ctrl->rlc_status[lcid].bytes_in_buffer == 0)
        continue;

      sched_ctrl->dl_pdus_total += sched_ctrl->rlc_status[lcid].pdus_in_buffer;
      sched_ctrl->num_total_bytes += sched_ctrl->rlc_status[lcid].bytes_in_buffer;
      LOG_D(MAC,
            "[gNB %d][%4d.%2d] %s%d->DLSCH, RLC status for UE %d: %d bytes in buffer, total DL buffer size = %d bytes, %d total PDU bytes, %s TA command\n",
            module_id,
            frame,
            slot,
            lcid < 4 ? "DCCH":"DTCH",
            lcid,
            UE->rnti,
            sched_ctrl->rlc_status[lcid].bytes_in_buffer,
            sched_ctrl->num_total_bytes,
            sched_ctrl->dl_pdus_total,
            sched_ctrl->ta_apply ? "send":"do not send");
    }
  }
}

typedef struct {
  int bwpStart;
  int bwpSize;
} dl_bwp_info_t;

static dl_bwp_info_t get_bwp_start_size(gNB_MAC_INST *mac, NR_UE_info_t *UE)
{
  NR_UE_DL_BWP_t *dl_bwp = &UE->current_DL_BWP;
  NR_UE_sched_ctrl_t *sched_ctrl = &UE->UE_sched_ctrl;
  // UE is scheduled in a set of contiguously allocated resource blocks within the active bandwidth part of size N_BWP PRBs
  // except for the case when DCI format 1_0 is decoded in any common search space
  // in which case the size of CORESET 0 shall be used if CORESET 0 is configured for the cell
  // and the size of initial DL bandwidth part shall be used if CORESET 0 is not configured for the cell.
  // TS 38.214 Section 5.1.2.2.2
  dl_bwp_info_t bwp_info;
  bwp_info.bwpSize = dl_bwp->BWPSize;
  bwp_info.bwpStart = dl_bwp->BWPStart;
  if (sched_ctrl->search_space->searchSpaceType->present == NR_SearchSpace__searchSpaceType_PR_common &&
      dl_bwp->dci_format == NR_DL_DCI_FORMAT_1_0) {
    if (mac->cset0_bwp_size != 0) {
      bwp_info.bwpStart = mac->cset0_bwp_start;
      bwp_info.bwpSize = mac->cset0_bwp_size;
    }
    else {
      // TODO this is not entirely correct
      // start would be the start of CORESET not of the initial BWP
      bwp_info.bwpStart = UE->sc_info.initial_dl_BWPStart;
      bwp_info.bwpSize = UE->sc_info.initial_dl_BWPSize;
    }
  }
  return bwp_info;
}


typedef struct UEsched_s {
  float coef;
  NR_UE_info_t * UE;
} UEsched_t;

void create_ue_sched_list( UEsched_t *UE_sched,
  mac_ctrl_msg_t mac_ctrl_msg,
  NR_UE_info_t **UE_list ) 
{
  for(int i=0; i<mac_ctrl_msg.num_users; ++i) {
    int cur_rnti = mac_ctrl_msg.resource_alloc[i].user_id;
    UE_iterator(UE_list, UE) {
      if(UE->rnti == cur_rnti) {
        UE_sched[i].UE = UE;
        break;
      }
    }
  }
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

static void do_resource_allocation(mac_ctrl_msg_t mac_ctrl_msg, module_id_t module_id,
                  frame_t frame,
                  slot_t slot,
                  NR_UE_info_t **UE_list,
                  int max_num_ue,
                  int num_beams,
                  int n_rb_sched[num_beams])
{
  gNB_MAC_INST *mac = RC.nrmac[module_id];
  NR_ServingCellConfigCommon_t *scc=mac->common_channels[0].ServingCellConfigCommon;
  // UEs that could be scheduled
  UEsched_t UE_sched[MAX_MOBILES_PER_GNB + 1] = {0};
  int remainUEs[num_beams];
  for (int i = 0; i < num_beams; i++)
    remainUEs[i] = max_num_ue;
  int curUE = 0;
  int CC_id = 0;
  int slots_per_frame = mac->frame_structure.numb_slots_frame;

  create_ue_sched_list(UE_sched, mac_ctrl_msg, UE_list);

  UEsched_t *iterator = UE_sched;

  const int min_rbSize = 5;

  
  while (iterator->UE != NULL) {

    NR_UE_sched_ctrl_t *sched_ctrl = &iterator->UE->UE_sched_ctrl;
    const uint16_t rnti = iterator->UE->rnti;

    NR_UE_DL_BWP_t *dl_bwp = &iterator->UE->current_DL_BWP;
    NR_UE_UL_BWP_t *ul_bwp = &iterator->UE->current_UL_BWP;

    if (sched_ctrl->available_dl_harq.head < 0) {
      LOG_D(NR_MAC, "[UE %04x][%4d.%2d] UE has no free DL HARQ process, skipping\n",
            iterator->UE->rnti,
            frame,
            slot);
      iterator++;
      continue;
    }

    NR_beam_alloc_t beam = beam_allocation_procedure(&mac->beam_info, frame, slot, iterator->UE->UE_beam_index, slots_per_frame);

    if (beam.idx < 0) {
      // no available beam
      iterator++;
      continue;
    }
    if (remainUEs[beam.idx] == 0 || n_rb_sched[beam.idx] < min_rbSize) {
      reset_beam_status(&mac->beam_info, frame, slot, iterator->UE->UE_beam_index, slots_per_frame, beam.new_beam);
      iterator++;
      continue;
    }

    NR_sched_pdsch_t *sched_pdsch = &sched_ctrl->sched_pdsch;
    sched_pdsch->dl_harq_pid = sched_ctrl->available_dl_harq.head;

    /* MCS has been set above */
    sched_pdsch->mcs = get_mcs_from_ctrl(rnti, mac_ctrl_msg);
    sched_pdsch->time_domain_allocation = get_dl_tda(mac, slot);
    AssertFatal(sched_pdsch->time_domain_allocation>=0,"Unable to find PDSCH time domain allocation in list\n");

    const int coresetid = sched_ctrl->coreset->controlResourceSetId;
    sched_pdsch->tda_info = get_dl_tda_info(dl_bwp,
                                            sched_ctrl->search_space->searchSpaceType->present,
                                            sched_pdsch->time_domain_allocation,
                                            scc->dmrs_TypeA_Position,
                                            1,
                                            TYPE_C_RNTI_,
                                            coresetid,
                                            false);
    AssertFatal(sched_pdsch->tda_info.valid_tda, "Invalid TDA from get_dl_tda_info\n");

    NR_tda_info_t *tda_info = &sched_pdsch->tda_info;

    const uint16_t slbitmap = SL_to_bitmap(tda_info->startSymbolIndex, tda_info->nrOfSymbols);

    uint16_t *rballoc_mask = mac->common_channels[CC_id].vrb_map[beam.idx];
    dl_bwp_info_t bwp_info = get_bwp_start_size(mac, iterator->UE);
    int rbStart = 0; // WRT BWP start
    int rbStop = bwp_info.bwpSize - 1;
    int bwp_start = bwp_info.bwpStart;
    // Freq-demain allocation
    while (rbStart < rbStop && (rballoc_mask[rbStart + bwp_start] & slbitmap))
      rbStart++;

    uint16_t max_rbSize = 1;

    while (rbStart + max_rbSize <= rbStop && !(rballoc_mask[rbStart + max_rbSize + bwp_start] & slbitmap))
      max_rbSize++;

    int ctrl_rbSize = get_rbSize_from_ctrl(rnti, mac_ctrl_msg);
    
    max_rbSize = (max_rbSize > ctrl_rbSize) ? max_rbSize : ctrl_rbSize;
    if (max_rbSize < min_rbSize) {
      LOG_D(NR_MAC,
            "(%d.%d) Cannot schedule RNTI %04x, rbStart %d, rbSize %d, rbStop %d\n",
            frame,
            slot,
            rnti,
            rbStart,
            max_rbSize,
            rbStop);
      iterator++;
      continue;
    }

    int CCEIndex = get_cce_index(mac,
                                 CC_id,
                                 slot,
                                 iterator->UE->rnti,
                                 &sched_ctrl->aggregation_level,
                                 beam.idx,
                                 sched_ctrl->search_space,
                                 sched_ctrl->coreset,
                                 &sched_ctrl->sched_pdcch,
                                 false,
                                 sched_ctrl->pdcch_cl_adjust);
    if (CCEIndex < 0) {
      LOG_D(NR_MAC, "[UE %04x][%4d.%2d] could not find free CCE for DL DCI\n", rnti, frame, slot);
      reset_beam_status(&mac->beam_info, frame, slot, iterator->UE->UE_beam_index, slots_per_frame, beam.new_beam);
      iterator++;
      continue;
    }

    /* Find PUCCH occasion: if it fails, undo CCE allocation (undoing PUCCH
    * allocation after CCE alloc fail would be more complex) */

    int alloc = -1;
    if (!get_FeedbackDisabled(iterator->UE->sc_info.downlinkHARQ_FeedbackDisabled_r17, sched_pdsch->dl_harq_pid)) {
      int r_pucch = nr_get_pucch_resource(sched_ctrl->coreset, ul_bwp->pucch_Config, CCEIndex);
      alloc = nr_acknack_scheduling(mac, iterator->UE, frame, slot, iterator->UE->UE_beam_index, r_pucch, 0);
      if (alloc < 0) {
        LOG_D(NR_MAC, "[UE %04x][%4d.%2d] could not find PUCCH for DL DCI\n", rnti, frame, slot);
        reset_beam_status(&mac->beam_info, frame, slot, iterator->UE->UE_beam_index, slots_per_frame, beam.new_beam);
        iterator++;
        continue;
      }
    }

    sched_ctrl->cce_index = CCEIndex;
    fill_pdcch_vrb_map(mac, CC_id, &sched_ctrl->sched_pdcch, CCEIndex, sched_ctrl->aggregation_level, beam.idx);

    sched_pdsch->dmrs_parms = get_dl_dmrs_params(scc, dl_bwp, tda_info, sched_pdsch->nrOfLayers);
    sched_pdsch->Qm = nr_get_Qm_dl(sched_pdsch->mcs, dl_bwp->mcsTableIdx);
    sched_pdsch->R = nr_get_code_rate_dl(sched_pdsch->mcs, dl_bwp->mcsTableIdx);
    sched_pdsch->pucch_allocation = alloc;
    uint32_t TBS = 0;
    uint16_t rbSize;
    // Fix me: currently, the RLC does not give us the total number of PDUs
    // awaiting. Therefore, for the time being, we put a fixed overhead of 12
    // (for 4 PDUs) and optionally + 2 for TA. Once RLC gives the number of
    // PDUs, we replace with 3 * numPDUs
    const int oh = 3 * 4 + 2 * (frame == (sched_ctrl->ta_frame + 100) % 1024);
    //const int oh = 3 * sched_ctrl->dl_pdus_total + 2 * (frame == (sched_ctrl->ta_frame + 100) % 1024);
    nr_find_nb_rb(sched_pdsch->Qm,
                  sched_pdsch->R,
                  1, // no transform precoding for DL
                  sched_pdsch->nrOfLayers,
                  tda_info->nrOfSymbols,
                  sched_pdsch->dmrs_parms.N_PRB_DMRS * sched_pdsch->dmrs_parms.N_DMRS_SLOT,
                  sched_ctrl->num_total_bytes + oh,
                  min_rbSize,
                  max_rbSize,
                  &TBS,
                  &rbSize);
    sched_pdsch->rbSize = rbSize;
    sched_pdsch->rbStart = rbStart;
    sched_pdsch->tb_size = TBS;
    /* transmissions: directly allocate */
    n_rb_sched[beam.idx] -= sched_pdsch->rbSize;

    for (int rb = bwp_start; rb < sched_pdsch->rbSize; rb++)
      rballoc_mask[rb + sched_pdsch->rbStart] |= slbitmap;

    remainUEs[beam.idx]--;
    iterator++;
  }
}

static void do_dlsch_preprocessing(module_id_t module_id, frame_t frame, slot_t slot, mac_ctrl_msg_t mac_ctrl_msg) {
  printf("CUSTOM (xApp driven) nr_dlsch_preprocessor: frame %d, slot %d\n",frame, slot);
  gNB_MAC_INST *mac = RC.nrmac[module_id];
  NR_UEs_t *UE_info = &mac->UE_info;

  if (UE_info->connected_ue_list[0] == NULL)
    return;

  NR_ServingCellConfigCommon_t *scc = mac->common_channels[0].ServingCellConfigCommon;
  int bw = scc->downlinkConfigCommon->frequencyInfoDL->scs_SpecificCarrierList.list.array[0]->carrierBandwidth;
  int num_beams = mac->beam_info.beam_allocation ? mac->beam_info.beams_per_period : 1;
  int n_rb_sched[num_beams];
  for (int i = 0; i < num_beams; i++)
    n_rb_sched[i] = bw;

  /* Retrieve amount of data to send for this UE */
  nr_store_dlsch_buffer(module_id, frame, slot);

  int average_agg_level = 4; // TODO find a better estimation
  int max_sched_ues = bw / (average_agg_level * NR_NB_REG_PER_CCE);

  // FAPI cannot handle more than MAX_DCI_CORESET DCIs
  max_sched_ues = min(max_sched_ues, MAX_DCI_CORESET);

  /* do per UE resource allocation according to MAC control message. */
  do_resource_allocation(
    mac_ctrl_msg,
    mod_id, frame, slot, UE_info->connected_ue_list, max_sched_ues, num_beams, n_rb_sched
  );
}

sm_ag_if_ans_t write_ctrl_mac_sm(void const* data)
{
  // This works. We actually do come here when E2 node gets a control request.
  //printf("\n\n\n\n\nIN THE CTRL CALLBACK!\n\n\n\n\n");
  assert(data != NULL);
  
  /*Get the preprocessor parameters ready*/
  frame_t frame = RC.ru[0]->proc.frame_tx; slot_t slot = RC.ru[0]->proc.tti_tx;
  module_id_t module_id = 0;

  /*Get the MAC Control message for the scheduling decisions*/
  mac_ctrl_req_data_t *mac_data = data;
  mac_ctrl_msg_t mac_ctrl_msg = mac_data->msg;

  /*Run the preprocessor for dlsch resource allocation decision.*/
  do_dlsch_preprocessing(module_id, frame, slot, mac_ctrl_msg);
  //sem_post(&custom_scheduler);

  sm_ag_if_ans_t ans = {0};
  return ans;
}

