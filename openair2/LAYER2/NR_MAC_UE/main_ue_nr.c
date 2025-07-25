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

/* \file main_ue_nr.c
 * \brief top init of Layer 2
 * \author R. Knopp, K.H. HSU
 * \date 2018
 * \version 0.1
 * \company Eurecom / NTUST
 * \email: knopp@eurecom.fr, kai-hsiang.hsu@eurecom.fr
 * \note
 * \warning
 */

//#include "defs.h"
#include "mac_proto.h"
#include "radio/COMMON/common_lib.h"
#include "assertions.h"
#include "executables/nr-uesoftmodem.h"
#include "nr_rlc/nr_rlc_oai_api.h"
#include "RRC/NR_UE/rrc_proto.h"
#include <pthread.h>
static NR_UE_MAC_INST_t *nr_ue_mac_inst; 

void send_srb0_rrc(int ue_id, const uint8_t *sdu, sdu_size_t sdu_len, void *data)
{
  AssertFatal(sdu_len > 0 && sdu_len < CCCH_SDU_SIZE, "invalid CCCH SDU size %d\n", sdu_len);

  MessageDef *message_p = itti_alloc_new_message(TASK_MAC_UE, 0, NR_RRC_MAC_CCCH_DATA_IND);
  memset(NR_RRC_MAC_CCCH_DATA_IND(message_p).sdu, 0, sdu_len);
  memcpy(NR_RRC_MAC_CCCH_DATA_IND(message_p).sdu, sdu, sdu_len);
  NR_RRC_MAC_CCCH_DATA_IND(message_p).sdu_size = sdu_len;
  itti_send_msg_to_task(TASK_RRC_NRUE, ue_id, message_p);
}

void nr_ue_init_mac(NR_UE_MAC_INST_t *mac)
{
  LOG_I(NR_MAC, "[UE%d] Initializing MAC\n", mac->ue_id);
  nr_ue_reset_sync_state(mac);
  mac->get_sib1 = false;
  for (int i = 0; i < MAX_SI_GROUPS; i++)
    mac->get_otherSI[i] = false;
  memset(&mac->phy_config, 0, sizeof(mac->phy_config));
  mac->si_SchedInfo.si_window_start = -1;
  mac->servCellIndex = 0;
  mac->harq_ACK_SpatialBundlingPUCCH = false;
  mac->harq_ACK_SpatialBundlingPUSCH = false;
  mac->uecap_maxMIMO_PDSCH_layers = 0;
  mac->uecap_maxMIMO_PUSCH_layers_cb = 0;
  mac->uecap_maxMIMO_PUSCH_layers_nocb = 0;
  mac->p_Max = INT_MIN;
  mac->p_Max_alt = INT_MIN;
  mac->msg3_C_RNTI = false;
  mac->ntn_ta.ntn_params_changed = false;
  reset_mac_inst(mac);

  // need to inizialize because might not been setup (optional timer)
  nr_timer_stop(&mac->scheduling_info.sr_DelayTimer);

  memset(&mac->ssb_measurements, 0, sizeof(mac->ssb_measurements));
  memset(&mac->ul_time_alignment, 0, sizeof(mac->ul_time_alignment));
  memset(&mac->ssb_list, 0, sizeof(mac->ssb_list));

  for (int i = 0; i < NR_MAX_SR_ID; i++)
    memset(&mac->scheduling_info.sr_info[i], 0, sizeof(mac->scheduling_info.sr_info[i]));

  mac->pucch_power_control_initialized = false;
  mac->pusch_power_control_initialized = false;
}

void nr_ue_mac_default_configs(NR_UE_MAC_INST_t *mac)
{
  // default values as defined in 38.331 sec 9.2.2

  // sf80 default for retxBSR_Timer sf10 for periodicBSR_Timer
  int mu = mac->current_UL_BWP ? mac->current_UL_BWP->scs : get_softmodem_params()->numerology;
  int subframes_per_slot = get_slots_per_frame_from_scs(mu) / 10;
  nr_timer_setup(&mac->scheduling_info.retxBSR_Timer, 80 * subframes_per_slot, 1); // 1 slot update rate
  nr_timer_setup(&mac->scheduling_info.periodicBSR_Timer, 10 * subframes_per_slot, 1); // 1 slot update rate

  mac->scheduling_info.phr_info.is_configured = true;
  mac->scheduling_info.phr_info.PathlossChange_db = 1;
  nr_timer_setup(&mac->scheduling_info.phr_info.periodicPHR_Timer, 10 * subframes_per_slot, 1);
  nr_timer_setup(&mac->scheduling_info.phr_info.prohibitPHR_Timer, 10 * subframes_per_slot, 1);
}

void nr_ue_send_synch_request(NR_UE_MAC_INST_t *mac, module_id_t module_id, int cc_id, const fapi_nr_synch_request_t *sync_req)
{
  // Sending to PHY a request to resync
  mac->synch_request.Mod_id = module_id;
  mac->synch_request.CC_id = cc_id;
  mac->synch_request.synch_req = *sync_req;
  mac->if_module->synch_request(&mac->synch_request);
}

void nr_ue_reset_sync_state(NR_UE_MAC_INST_t *mac)
{
  // reset synchornization status
  mac->state = UE_NOT_SYNC;
  mac->ra.ra_state = nrRA_UE_IDLE;
}

NR_UE_L2_STATE_t nr_ue_get_sync_state(module_id_t mod_id)
{
  NR_UE_MAC_INST_t *mac = get_mac_inst(mod_id);
  return mac->state;
}

NR_UE_MAC_INST_t *nr_l2_init_ue(int nb_inst)
{
  //init mac here
  nr_ue_mac_inst = (NR_UE_MAC_INST_t *)calloc(nb_inst, sizeof(NR_UE_MAC_INST_t));
  AssertFatal(nr_ue_mac_inst, "Couldn't allocate %d instances of MAC module\n", nb_inst);

  for (int j = 0; j < nb_inst; j++) {
    NR_UE_MAC_INST_t *mac = &nr_ue_mac_inst[j];
    mac->ue_id = j;
    nr_ue_init_mac(mac);
    int ret = pthread_mutex_init(&mac->if_mutex, NULL);
    AssertFatal(ret == 0, "Mutex init failed\n");
    nr_ue_mac_default_configs(mac);
    if (IS_SA_MODE(get_softmodem_params()))
      ue_init_config_request(mac, get_slots_per_frame_from_scs(get_softmodem_params()->numerology));
  }

  int rc = nr_rlc_module_init(NR_RLC_OP_MODE_UE);
  AssertFatal(rc == 0, "Could not initialize RLC layer\n");

  for (int j = 0; j < nb_inst; j++) {
    nr_rlc_activate_srb0(j, NULL, send_srb0_rrc);
  }

  return (nr_ue_mac_inst);
}

NR_UE_MAC_INST_t *get_mac_inst(module_id_t module_id)
{
  NR_UE_MAC_INST_t *mac = &nr_ue_mac_inst[(int)module_id];
  AssertFatal(mac, "Couldn't get MAC inst %d\n", module_id);
  AssertFatal(mac->ue_id == module_id, "MAC ID %d doesn't match with input %d\n", mac->ue_id, module_id);
  return mac;
}

void reset_mac_inst(NR_UE_MAC_INST_t *nr_mac)
{
  // MAC reset according to 38.321 Section 5.12

  // initialize Bj for each logical channel to zero
  // TODO reset also other status variables of LC, is this ok?
  for (int i = 0; i < NR_MAX_NUM_LCID; i++) {
    LOG_D(NR_MAC, "Applying default logical channel config for LCID %d\n", i);
    nr_mac->scheduling_info.lc_sched_info[i].Bj = 0;
    nr_mac->scheduling_info.lc_sched_info[i].LCID_buffer_remain = 0;
  }

  // TODO stop all running timers
  for (int i = 0; i < NR_MAX_NUM_LCID; i++) {
    nr_mac->scheduling_info.lc_sched_info[i].Bj = 0;
    nr_timer_stop(&nr_mac->scheduling_info.lc_sched_info[i].Bj_timer);
  }
  if (nr_mac->data_inactivity_timer)
    nr_timer_stop(nr_mac->data_inactivity_timer);
  nr_timer_stop(&nr_mac->time_alignment_timer);
  nr_timer_stop(&nr_mac->scheduling_info.sr_DelayTimer);
  nr_timer_stop(&nr_mac->scheduling_info.retxBSR_Timer);
  nr_timer_stop(&nr_mac->ra.response_window_timer);
  nr_timer_stop(&nr_mac->ra.RA_backoff_timer);
  nr_timer_stop(&nr_mac->ra.contention_resolution_timer);
  for (int i = 0; i < NR_MAX_SR_ID; i++)
    nr_timer_stop(&nr_mac->scheduling_info.sr_info[i].prohibitTimer);

  // consider all timeAlignmentTimers as expired and perform the corresponding actions in clause 5.2
  handle_time_alignment_timer_expired(nr_mac);

  // set the NDIs for all uplink HARQ processes to the value 0
  for (int k = 0; k < NR_MAX_HARQ_PROCESSES; k++)
    nr_mac->ul_harq_info[k].last_ndi = -1; // initialize to invalid value

  // stop any ongoing RACH procedure
  if (nr_mac->ra.RA_active) {
    nr_mac->ra.ra_state = nrRA_UE_IDLE;
    nr_mac->ra.RA_active = false;
  }

  // discard explicitly signalled contention-free Random Access Resources
  // TODO not sure what needs to be done here

  // flush Msg3 buffer
  free_and_zero(nr_mac->ra.Msg3_buffer);

  // cancel any triggered Scheduling Request procedure
  for (int i = 0; i < NR_MAX_SR_ID; i++) {
    nr_mac->scheduling_info.sr_info[i].pending = false;
    nr_mac->scheduling_info.sr_info[i].counter = 0;
  }

  // cancel any triggered Buffer Status Reporting procedure
  nr_mac->scheduling_info.BSR_reporting_active = NR_BSR_TRIGGER_NONE;

  // cancel any triggered Power Headroom Reporting procedure
  nr_mac->scheduling_info.phr_info.phr_reporting = 0;
  nr_mac->scheduling_info.phr_info.was_mac_reset = true;

  // flush the soft buffers for all DL HARQ processes
  memset(nr_mac->dl_harq_info, 0, sizeof(nr_mac->dl_harq_info));

  // for each DL HARQ process, consider the next received transmission for a TB as the very first transmission
  for (int k = 0; k < NR_MAX_HARQ_PROCESSES; k++)
    nr_mac->dl_harq_info[k].last_ndi = -1; // initialize to invalid value

  // release, if any, Temporary C-RNTI
  nr_mac->ra.t_crnti = 0;

  // reset BFI_COUNTER
  // TODO beam failure procedure not implemented
}

void release_mac_configuration(NR_UE_MAC_INST_t *mac, NR_UE_MAC_reset_cause_t cause)
{
  NR_UE_ServingCell_Info_t *sc = &mac->sc_info;
  // if cause is Re-establishment, release spCellConfig only
  if (cause == GO_TO_IDLE) {
    asn1cFreeStruc(asn_DEF_NR_MIB, mac->mib);
    asn1cFreeStruc(asn_DEF_NR_SearchSpace, mac->search_space_zero);
    asn1cFreeStruc(asn_DEF_NR_ControlResourceSet, mac->coreset0);
    asn_sequence_empty(&mac->si_SchedInfo.si_SchedInfo_list);
    asn1cFreeStruc(asn_DEF_NR_TDD_UL_DL_ConfigCommon, mac->tdd_UL_DL_ConfigurationCommon);
    for (int i = mac->lc_ordered_list.count; i > 0 ; i--)
      asn_sequence_del(&mac->lc_ordered_list, i - 1, 1);
  }

  asn1cFreeStruc(asn_DEF_NR_CrossCarrierSchedulingConfig, sc->crossCarrierSchedulingConfig);
  asn1cFreeStruc(asn_DEF_NR_SRS_CarrierSwitching, sc->carrierSwitching);
  asn1cFreeStruc(asn_DEF_NR_UplinkConfig, sc->supplementaryUplink);
  asn1cFreeStruc(asn_DEF_NR_PDSCH_CodeBlockGroupTransmission, sc->pdsch_CGB_Transmission);
  asn1cFreeStruc(asn_DEF_NR_PUSCH_CodeBlockGroupTransmission, sc->pusch_CGB_Transmission);
  asn1cFreeStruc(asn_DEF_NR_CSI_MeasConfig, sc->csi_MeasConfig);
  asn1cFreeStruc(asn_DEF_NR_CSI_AperiodicTriggerStateList, sc->aperiodicTriggerStateList);
  asn1cFreeStruc(asn_DEF_NR_NTN_Config_r17, sc->ntn_Config_r17);
  asn1cFreeStruc(asn_DEF_NR_DownlinkHARQ_FeedbackDisabled_r17, sc->downlinkHARQ_FeedbackDisabled_r17);
  free(sc->xOverhead_PDSCH);
  free(sc->nrofHARQ_ProcessesForPDSCH);
  free(sc->nrofHARQ_ProcessesForPDSCH_v1700);
  free(sc->nrofHARQ_ProcessesForPUSCH_r17);
  free(sc->rateMatching_PUSCH);
  free(sc->xOverhead_PUSCH);
  free(sc->maxMIMO_Layers_PDSCH);
  free(sc->maxMIMO_Layers_PUSCH);
  memset(&mac->sc_info, 0, sizeof(mac->sc_info));

  mac->current_DL_BWP = NULL;
  mac->current_UL_BWP = NULL;

  // in case of re-establishment we don't need to release initial BWP config common
  int first_bwp_rel = 0; // first BWP to release
  if (cause == RE_ESTABLISHMENT) {
    first_bwp_rel = 1;
    // release dedicated BWP0 config
    NR_UE_DL_BWP_t *bwp = mac->dl_BWPs.array[0];
    NR_BWP_PDCCH_t *pdcch = &mac->config_BWP_PDCCH[0];
    for (int i = pdcch->list_Coreset.count; i > 0 ; i--)
      asn_sequence_del(&pdcch->list_Coreset, i - 1, 1);
    for (int i = pdcch->list_SS.count; i > 0 ; i--)
      asn_sequence_del(&pdcch->list_SS, i - 1, 1);
    asn1cFreeStruc(asn_DEF_NR_PDSCH_Config, bwp->pdsch_Config);
    NR_UE_UL_BWP_t *ubwp = mac->ul_BWPs.array[0];
    asn1cFreeStruc(asn_DEF_NR_PUCCH_Config, ubwp->pucch_Config);
    asn1cFreeStruc(asn_DEF_NR_SRS_Config, ubwp->srs_Config);
    asn1cFreeStruc(asn_DEF_NR_PUSCH_Config, ubwp->pusch_Config);
    mac->current_DL_BWP = bwp;
    mac->current_UL_BWP = ubwp;
    mac->sc_info.initial_dl_BWPSize = bwp->BWPSize;
    mac->sc_info.initial_dl_BWPStart = bwp->BWPStart;
    mac->sc_info.initial_ul_BWPSize = ubwp->BWPSize;
    mac->sc_info.initial_ul_BWPStart = ubwp->BWPStart;
  }

  for (int i = first_bwp_rel; i < mac->dl_BWPs.count; i++)
    release_dl_BWP(mac, i);
  for (int i = first_bwp_rel; i < mac->ul_BWPs.count; i++)
    release_ul_BWP(mac, i);

  memset(&mac->ssb_measurements, 0, sizeof(mac->ssb_measurements));
  memset(&mac->csirs_measurements, 0, sizeof(mac->csirs_measurements));
  memset(&mac->ul_time_alignment, 0, sizeof(mac->ul_time_alignment));
  for (int i = mac->TAG_list.count; i > 0 ; i--)
    asn_sequence_del(&mac->TAG_list, i - 1, 1);
}
