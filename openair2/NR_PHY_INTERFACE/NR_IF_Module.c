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

/*! \file openair2/NR_PHY_INTERFACE/NR_IF_Module.c
* \brief data structures for PHY/MAC interface modules
* \author EURECOM/NTUST
* \date 2018
* \version 0.1
* \company Eurecom, NTUST
* \email: raymond.knopp@eurecom.fr, kroempa@gmail.com
* \note
* \warning
*/

#include "openair2/NR_PHY_INTERFACE/NR_IF_Module.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "LAYER2/NR_MAC_gNB/mac_proto.h"
#include "common/platform_constants.h"
#include "common/utils/T/T.h"
#include "executables/softmodem-common.h"
#include "nfapi/oai_integration/gnb_ind_vars.h"
#include "nfapi/oai_integration/vendor_ext.h"
#include "nfapi_interface.h"
#include "openair2/NR_PHY_INTERFACE/nr_sched_response.h"
#include "openair2/PHY_INTERFACE/queue_t.h"
#include "utils.h"
#include "nfapi/oai_integration/nfapi_pnf.h"

#define MAX_IF_MODULES 100

static NR_IF_Module_t *nr_if_inst[MAX_IF_MODULES];

static void handle_nr_rach(NR_UL_IND_t *UL_info)
{
  if (NFAPI_MODE == NFAPI_MODE_PNF) {
    if (UL_info->rach_ind.number_of_pdus > 0) {
      LOG_D(PHY,"UL_info->UL_info->rach_ind.number_of_pdus:%d SFN/Slot:%d.%d \n", UL_info->rach_ind.number_of_pdus, UL_info->rach_ind.sfn,UL_info->rach_ind.slot);
      oai_nfapi_nr_rach_indication(&UL_info->rach_ind);
      UL_info->rach_ind.number_of_pdus = 0;
    }
    return;
  }
  if (UL_info->rach_ind.number_of_pdus) {
    LOG_D(MAC,
          "UL_info[Frame %d, Slot %d] Calling initiate_ra_proc RACH:SFN/SLOT:%d/%d\n",
          UL_info->frame,
          UL_info->slot,
          UL_info->rach_ind.sfn,
          UL_info->rach_ind.slot);
    for (int i = 0; i < UL_info->rach_ind.number_of_pdus; i++) {
      nfapi_nr_prach_indication_pdu_t *rach = UL_info->rach_ind.pdu_list + i;
      if (rach->num_preamble > 1) {
        LOG_E(MAC, "Not more than 1 preamble per RACH PDU supported, ignoring the rest\n");
      }
      nr_initiate_ra_proc(UL_info->module_id,
                          UL_info->CC_id,
                          UL_info->rach_ind.sfn,
                          UL_info->rach_ind.slot,
                          rach->preamble_list[0].preamble_index,
                          rach->freq_index,
                          rach->symbol_index,
                          rach->preamble_list[0].timing_advance,
                          rach->preamble_list[0].preamble_pwr);
    }
  }
}

static void handle_nr_uci(NR_UL_IND_t *UL_info)
{
  if(NFAPI_MODE == NFAPI_MODE_PNF) {
    if (UL_info->uci_ind.num_ucis > 0) {
      LOG_D(PHY,"PNF Sending UL_info->num_ucis:%d PDU_type: %d, SFN/SF:%d.%d \n", UL_info->uci_ind.num_ucis, UL_info->uci_ind.uci_list[0].pdu_type ,UL_info->frame, UL_info->slot);
      oai_nfapi_nr_uci_indication(&UL_info->uci_ind);
      UL_info->uci_ind.num_ucis = 0;
    }
    return;
  }

  const module_id_t mod_id = UL_info->module_id;
  const frame_t frame = UL_info->uci_ind.sfn;
  const slot_t slot = UL_info->uci_ind.slot;
  int num_ucis = UL_info->uci_ind.num_ucis;
  nfapi_nr_uci_t *uci_list = UL_info->uci_ind.uci_list;

  for (int i = 0; i < num_ucis; i++) {
    switch (uci_list[i].pdu_type) {
      case NFAPI_NR_UCI_PUSCH_PDU_TYPE:
        LOG_E(MAC, "%s(): unhandled NFAPI_NR_UCI_PUSCH_PDU_TYPE\n", __func__);
        break;

      case NFAPI_NR_UCI_FORMAT_0_1_PDU_TYPE: {
        const nfapi_nr_uci_pucch_pdu_format_0_1_t *uci_pdu = &uci_list[i].pucch_pdu_format_0_1;
        LOG_D(NR_MAC, "The received uci has sfn slot %d %d, num_ucis %d and pdu_size %d\n",
                UL_info->uci_ind.sfn, UL_info->uci_ind.slot, num_ucis, uci_list[i].pdu_size);
        handle_nr_uci_pucch_0_1(mod_id, frame, slot, uci_pdu);
        break;
      }

        case NFAPI_NR_UCI_FORMAT_2_3_4_PDU_TYPE: {
          const nfapi_nr_uci_pucch_pdu_format_2_3_4_t *uci_pdu = &uci_list[i].pucch_pdu_format_2_3_4;
          handle_nr_uci_pucch_2_3_4(mod_id, frame, slot, uci_pdu);
          break;
        }
      LOG_D(MAC, "UCI handled \n");
    }
  }

  UL_info->uci_ind.num_ucis = 0;

}

struct sfn_slot {
  int sfn;
  int slot;
};
static bool crc_sfn_slot_matcher(void *wanted, void *candidate)
{
  nfapi_nr_p7_message_header_t *msg = candidate;
  struct sfn_slot *sfn_sf = (struct sfn_slot *)wanted;

  switch (msg->message_id) {
    case NFAPI_NR_PHY_MSG_TYPE_CRC_INDICATION: {
      nfapi_nr_crc_indication_t *ind = candidate;
      return sfn_sf->sfn == ind->sfn && sfn_sf->slot == ind->slot;
    }

    default:
      LOG_E(NR_MAC, "sfn_slot_match bad ID: %d\n", msg->message_id);
  }
  return false;
}

static void handle_nr_ulsch(NR_UL_IND_t *UL_info)
{
  if(NFAPI_MODE == NFAPI_MODE_PNF) {
    if (UL_info->crc_ind.number_crcs > 0) {
      LOG_D(PHY,"UL_info->UL_info->crc_ind.number_crcs:%d CRC_IND:SFN/Slot:%d.%d\n", UL_info->crc_ind.number_crcs, UL_info->crc_ind.sfn, UL_info->crc_ind.slot);
      oai_nfapi_nr_crc_indication(&UL_info->crc_ind);
      UL_info->crc_ind.number_crcs = 0;
    }

    if (UL_info->rx_ind.number_of_pdus > 0) {
      LOG_D(PHY,"UL_info->rx_ind.number_of_pdus:%d RX_IND:SFN/Slot:%d.%d \n", UL_info->rx_ind.number_of_pdus, UL_info->rx_ind.sfn, UL_info->rx_ind.slot);
      oai_nfapi_nr_rx_data_indication(&UL_info->rx_ind);
      UL_info->rx_ind.number_of_pdus = 0;
    }
    return;
  }

  if (UL_info->rx_ind.number_of_pdus > 0 && UL_info->crc_ind.number_crcs > 0) {
    AssertFatal(UL_info->rx_ind.number_of_pdus == UL_info->crc_ind.number_crcs,
                "number_of_pdus %d, number_crcs %d\n",
                UL_info->rx_ind.number_of_pdus, UL_info->crc_ind.number_crcs);
    for (int i = 0; i < UL_info->rx_ind.number_of_pdus; i++) {
      const nfapi_nr_rx_data_pdu_t *rx = &UL_info->rx_ind.pdu_list[i];
      const nfapi_nr_crc_t *crc = &UL_info->crc_ind.crc_list[i];
      LOG_D(NR_PHY, "UL_info->crc_ind.pdu_list[%d].rnti:%04x "
                    "UL_info->rx_ind.pdu_list[%d].rnti:%04x\n",
                    i, crc->rnti, i, rx->rnti);
      AssertFatal(crc->rnti == rx->rnti, "mis-match between CRC RNTI %04x and RX RNTI %04x\n", crc->rnti, rx->rnti);

      AssertFatal(crc->harq_id == rx->harq_id,
                  "mis-match between CRC HARQ ID %04x and RX HARQ ID %04x\n",
                  crc->harq_id,
                  rx->harq_id);

      LOG_D(NR_MAC,
            "%4d.%2d Calling rx_sdu (CRC %s/tb_crc_status %d)\n",
            UL_info->frame,
            UL_info->slot,
            crc->tb_crc_status ? "error" : "ok",
            crc->tb_crc_status);

      /* if CRC passes, pass PDU, otherwise pass NULL as error indication */
      nr_rx_sdu(UL_info->module_id,
                UL_info->CC_id,
                UL_info->rx_ind.sfn,
                UL_info->rx_ind.slot,
                crc->rnti,
                crc->tb_crc_status ? NULL : rx->pdu,
                rx->pdu_length,
                crc->harq_id,
                crc->timing_advance,
                crc->ul_cqi,
                crc->rssi);
    }
    UL_info->rx_ind.number_of_pdus = 0;
    UL_info->crc_ind.number_crcs = 0;
  }
}

static void handle_nr_srs(NR_UL_IND_t *UL_info)
{
  if(NFAPI_MODE == NFAPI_MODE_PNF) {
    if (UL_info->srs_ind.number_of_pdus > 0) {
      LOG_D(PHY,"PNF Sending UL_info->srs_ind.number_of_pdus: %d, SFN/SF:%d.%d \n",
            UL_info->srs_ind.number_of_pdus, UL_info->frame, UL_info->slot);
      oai_nfapi_nr_srs_indication(&UL_info->srs_ind);
      UL_info->srs_ind.number_of_pdus = 0;
    }
    return;
  }

  const module_id_t module_id = UL_info->module_id;
  const frame_t frame = UL_info->srs_ind.sfn;
  const slot_t slot = UL_info->srs_ind.slot;
  const int num_srs = UL_info->srs_ind.number_of_pdus;
  nfapi_nr_srs_indication_pdu_t *srs_list = UL_info->srs_ind.pdu_list;

  // from here

  for (int i = 0; i < num_srs; i++) {
    nfapi_nr_srs_indication_pdu_t *srs_ind = &srs_list[i];
    LOG_D(NR_PHY, "(%d.%d) UL_info->srs_ind.pdu_list[%d].rnti: 0x%04x\n", frame, slot, i, srs_ind->rnti);
    handle_nr_srs_measurements(module_id,
                               frame,
                               slot,
                               srs_ind);
  }

  UL_info->srs_ind.number_of_pdus = 0;
}

static void free_unqueued_nfapi_indications(nfapi_nr_rach_indication_t *rach_ind,
                                            nfapi_nr_uci_indication_t *uci_ind,
                                            nfapi_nr_rx_data_indication_t *rx_ind,
                                            nfapi_nr_crc_indication_t *crc_ind) {
  if (rach_ind && rach_ind->number_of_pdus > 0)
  {
    free_and_zero(rach_ind->pdu_list);
    free_and_zero(rach_ind);
  }
  if (uci_ind && uci_ind->num_ucis > 0)
  {
    /* PUCCH fields (HARQ, SR) are freed in handle_nr_uci_pucch_0_1() and
     * handle_nr_uci_pucch_2_3_4() */
    free_and_zero(uci_ind->uci_list);
    free_and_zero(uci_ind);
  }
  if (rx_ind && rx_ind->number_of_pdus > 0)
  {
    for (int i = 0; i < rx_ind->number_of_pdus; ++i) {
      free_and_zero(rx_ind->pdu_list[i].pdu);
    }
    free_and_zero(rx_ind->pdu_list);
    free_and_zero(rx_ind);
  }
  if (crc_ind && crc_ind->number_crcs > 0)
  {
    free_and_zero(crc_ind->crc_list);
    free_and_zero(crc_ind);
  }
}

static void remove_crc_pdu(nfapi_nr_crc_indication_t *crc_ind, int index) {
  AssertFatal(index >= 0, "Invalid index %d\n", index);
  AssertFatal(index < crc_ind->number_crcs, "Invalid index %d\n", index);
  AssertFatal(crc_ind->number_crcs > 0, "Invalid crc_ind->number_crcs %d\n", crc_ind->number_crcs);

  memmove(crc_ind->crc_list + index,
          crc_ind->crc_list + index + 1,
          sizeof(*crc_ind->crc_list) * (crc_ind->number_crcs - index - 1));
  crc_ind->number_crcs--;
}

static void remove_rx_pdu(nfapi_nr_rx_data_indication_t *rx_ind, int index) {
  AssertFatal(index >= 0, "Invalid index %d\n", index);
  AssertFatal(index < rx_ind->number_of_pdus, "Invalid index %d\n", index);
  AssertFatal(rx_ind->number_of_pdus > 0, "Invalid rx_ind->number_of_pdus %d\n", rx_ind->number_of_pdus);

  memmove(rx_ind->pdu_list + index,
          rx_ind->pdu_list + index + 1,
          sizeof(*rx_ind->pdu_list) * (rx_ind->number_of_pdus - index - 1));
  rx_ind->number_of_pdus--;
}

static bool crc_ind_has_rnti(nfapi_nr_crc_indication_t *crc_ind, uint16_t rnti) {
  for (int i = 0; i < crc_ind->number_crcs; i++) {
    if (rnti == crc_ind->crc_list[i].rnti) {
      return true;
    }
  }
  return false;
}

static bool rx_ind_has_rnti(nfapi_nr_rx_data_indication_t *rx_ind, uint16_t rnti) {
  for (int i = 0; i < rx_ind->number_of_pdus; i++) {
    if (rnti == rx_ind->pdu_list[i].rnti) {
      return true;
    }
  }
  return false;
}

static void match_crc_rx_pdu(nfapi_nr_rx_data_indication_t *rx_ind, nfapi_nr_crc_indication_t *crc_ind) {
  AssertFatal(crc_ind->number_crcs > 0 &&  rx_ind->number_of_pdus > 0,
              "Invalid number of crc_ind->number_crcs %d or rx_ind->number_of_pdus %d\n",
              crc_ind->number_crcs, rx_ind->number_of_pdus);
  if (crc_ind->number_crcs > rx_ind->number_of_pdus) {
    int num_unmatched_crcs = 0;
    nfapi_nr_crc_indication_t *crc_ind_unmatched = calloc(1, sizeof(*crc_ind_unmatched));
    crc_ind_unmatched->header = crc_ind->header;
    crc_ind_unmatched->sfn = crc_ind->sfn;
    crc_ind_unmatched->slot = crc_ind->slot;
    crc_ind_unmatched->number_crcs = crc_ind->number_crcs - rx_ind->number_of_pdus;
    crc_ind_unmatched->crc_list = calloc(crc_ind_unmatched->number_crcs, sizeof(nfapi_nr_crc_t));
    for (int i = 0; i < crc_ind->number_crcs; i++) {
      if (!rx_ind_has_rnti(rx_ind, crc_ind->crc_list[i].rnti)) {
          LOG_I(NR_MAC, "crc_ind->crc_list[%d].rnti %x does not match any rx_ind pdu rnti\n",
                i, crc_ind->crc_list[i].rnti);
          crc_ind_unmatched->crc_list[num_unmatched_crcs] = crc_ind->crc_list[i];
          num_unmatched_crcs++;
          remove_crc_pdu(crc_ind, i);
      }
      if (crc_ind->number_crcs == rx_ind->number_of_pdus) {
        break;
      }
    }
    AssertFatal(crc_ind_unmatched->number_crcs == num_unmatched_crcs, "crc_ind num_pdus %d doesnt match %d\n",
                crc_ind_unmatched->number_crcs, num_unmatched_crcs);
    if (!requeue(&gnb_crc_ind_queue, crc_ind_unmatched))
    {
      LOG_E(NR_PHY, "requeue failed for crc_ind_unmatched.\n");
      free_and_zero(crc_ind_unmatched->crc_list);
      free_and_zero(crc_ind_unmatched);
    }
  }
  else if (crc_ind->number_crcs < rx_ind->number_of_pdus) {
    int num_unmatched_rxs = 0;
    nfapi_nr_rx_data_indication_t *rx_ind_unmatched = calloc(1, sizeof(*rx_ind_unmatched));
    rx_ind_unmatched->header = rx_ind->header;
    rx_ind_unmatched->sfn = rx_ind->sfn;
    rx_ind_unmatched->slot = rx_ind->slot;
    rx_ind_unmatched->number_of_pdus = rx_ind->number_of_pdus - crc_ind->number_crcs;
    rx_ind_unmatched->pdu_list = calloc(rx_ind_unmatched->number_of_pdus, sizeof(nfapi_nr_pdu_t));
    for (int i = 0; i < rx_ind->number_of_pdus; i++) {
      if (!crc_ind_has_rnti(crc_ind, rx_ind->pdu_list[i].rnti)) {
        LOG_I(NR_MAC, "rx_ind->pdu_list[%d].rnti %x does not match any crc_ind pdu rnti\n", i, rx_ind->pdu_list[i].rnti);
        rx_ind_unmatched->pdu_list[num_unmatched_rxs] = rx_ind->pdu_list[i];
        num_unmatched_rxs++;
        remove_rx_pdu(rx_ind, i);
      }
      if (rx_ind->number_of_pdus == crc_ind->number_crcs) {
        break;
      }
    }
    AssertFatal(rx_ind_unmatched->number_of_pdus == num_unmatched_rxs, "rx_ind num_pdus %d doesnt match %d\n",
                rx_ind_unmatched->number_of_pdus, num_unmatched_rxs);
    if (!requeue(&gnb_rx_ind_queue, rx_ind_unmatched))
    {
      LOG_E(NR_PHY, "requeue failed for rx_ind_unmatched.\n");
      free_and_zero(rx_ind_unmatched->pdu_list);
      free_and_zero(rx_ind_unmatched);
    }
  }
  else {
    LOG_E(NR_MAC, "The number of crc pdus %d = the number of rx pdus %d\n",
          crc_ind->number_crcs, rx_ind->number_of_pdus);
  }
}

extern void handle_nr_slot_ind(uint16_t sfn, uint16_t slot);
static void pnf_send_slot_ind(module_id_t module_id, int CC_id, int frame, int slot)
{
  (void)module_id;
  (void)CC_id;
  handle_nr_slot_ind(frame, slot);
}

static void run_scheduler_monolithic(module_id_t module_id, int CC_id, int frame, int slot)
{

  // printf("\n\n\n\nrun_scheduler_monolithic(): HERE!!");
  // printf("\tframe=%d slot=%d\n\n\n\n", frame, slot);
  NR_IF_Module_t *ifi = nr_if_inst[module_id];

  LOG_D(NR_MAC, "Calling scheduler for %d.%d\n", frame, slot);
  NR_Sched_Rsp_t *sched_info = allocate_sched_response();

  // clear UL DCI prior to handling ULSCH
  sched_info->UL_dci_req.numPdus = 0;
  gNB_dlsch_ulsch_scheduler(module_id, frame, slot, sched_info);
  ifi->CC_mask = 0;
  sched_info->module_id = module_id;
  sched_info->CC_id = CC_id;
  sched_info->frame = frame;
  sched_info->slot = slot;
  /*
  sched_info->DL_req      = &mac->DL_req[CC_id];
  sched_info->UL_dci_req  = &mac->UL_dci_req[CC_id];

  sched_info->UL_tti_req  = mac->UL_tti_req[CC_id];

  sched_info->TX_req      = &mac->TX_req[CC_id];
  */
#ifdef DUMP_FAPI
  dump_dl(sched_info);
#endif

  AssertFatal(ifi->NR_Schedule_response != NULL, "nr_schedule_response is null (mod %d, cc %d)\n", module_id, CC_id);
  ifi->NR_Schedule_response(sched_info);

  LOG_D(NR_PHY,
        "NR_Schedule_response: SFN SLOT:%d %d dl_pdus:%d\n",
        sched_info->frame,
        sched_info->slot,
        sched_info->DL_req.dl_tti_request_body.nPDUs);
}

void NR_UL_indication(NR_UL_IND_t *UL_info) {
  AssertFatal(UL_info!=NULL,"UL_info is null\n");
  module_id_t module_id = UL_info->module_id;
  int CC_id = UL_info->CC_id;

  LOG_D(NR_PHY,"SFN/SLOT:%d.%d module_id:%d CC_id:%d UL_info[rach_pdus:%zu rx_ind:%zu crcs:%zu]\n",
        UL_info->frame, UL_info->slot,
        module_id, CC_id,
        gnb_rach_ind_queue.num_items,
        gnb_rx_ind_queue.num_items,
        gnb_crc_ind_queue.num_items);

  nfapi_nr_rach_indication_t *rach_ind = NULL;
  nfapi_nr_uci_indication_t *uci_ind = NULL;
  nfapi_nr_rx_data_indication_t *rx_ind = NULL;
  nfapi_nr_crc_indication_t *crc_ind = NULL;
  if (NFAPI_MODE == NFAPI_MODE_VNF || NFAPI_MODE == NFAPI_MODE_AERIAL)
  {
    if (gnb_rach_ind_queue.num_items > 0) {
      LOG_D(NR_MAC, "gnb_rach_ind_queue size = %zu\n", gnb_rach_ind_queue.num_items);
      rach_ind = get_queue(&gnb_rach_ind_queue);
      AssertFatal(rach_ind->number_of_pdus > 0, "Invalid number of PDUs\n");
      UL_info->rach_ind = *rach_ind;
    }
    if (gnb_uci_ind_queue.num_items > 0) {
      LOG_D(NR_MAC, "gnb_uci_ind_queue size = %zu\n", gnb_uci_ind_queue.num_items);
      uci_ind = get_queue(&gnb_uci_ind_queue);
      AssertFatal(uci_ind->num_ucis > 0, "Invalid number of PDUs\n");
      UL_info->uci_ind = *uci_ind;
    }
    if (gnb_rx_ind_queue.num_items > 0 && gnb_crc_ind_queue.num_items > 0) {
      LOG_D(NR_MAC, "gnb_rx_ind_queue size = %zu and gnb_crc_ind_queue size = %zu\n",
            gnb_rx_ind_queue.num_items, gnb_crc_ind_queue.num_items);
      rx_ind = get_queue(&gnb_rx_ind_queue);
      struct sfn_slot sfn_slot = {.sfn = rx_ind->sfn, .slot = rx_ind->slot};
      crc_ind = unqueue_matching(&gnb_crc_ind_queue,
                                 MAX_QUEUE_SIZE,
                                 crc_sfn_slot_matcher,
                                 &sfn_slot);
      if (!crc_ind) {
        LOG_I(NR_PHY, "No crc indication with the same SFN SLOT of rx indication %u %u\n", rx_ind->sfn, rx_ind->slot);
        requeue(&gnb_rx_ind_queue, rx_ind);
      }
      else {
        AssertFatal(rx_ind->number_of_pdus > 0, "Invalid number of PDUs\n");
        AssertFatal(crc_ind->number_crcs > 0, "Invalid number of PDUs\n");
        if (crc_ind->number_crcs != rx_ind->number_of_pdus)
          match_crc_rx_pdu(rx_ind, crc_ind);
        UL_info->rx_ind = *rx_ind;
        UL_info->crc_ind = *crc_ind;
      }
    }
  }

  if (UL_info->rach_ind.number_of_pdus > 0)
    handle_nr_rach(UL_info);
  handle_nr_uci(UL_info);
  handle_nr_ulsch(UL_info);
  handle_nr_srs(UL_info);

  if (NFAPI_MODE == NFAPI_MODE_VNF || NFAPI_MODE == NFAPI_MODE_AERIAL) {
    free_unqueued_nfapi_indications(rach_ind, uci_ind, rx_ind, crc_ind);
  }
}

NR_IF_Module_t *NR_IF_Module_init(int Mod_id) {
  AssertFatal(Mod_id<MAX_MODULES,"Asking for Module %d > %d\n",Mod_id,MAX_IF_MODULES);
  LOG_D(PHY, "Installing callbacks for IF_Module - UL_indication\n");

  if (nr_if_inst[Mod_id]==NULL) {
    nr_if_inst[Mod_id] = (NR_IF_Module_t*)malloc(sizeof(NR_IF_Module_t));
    memset((void*)nr_if_inst[Mod_id],0,sizeof(NR_IF_Module_t));

    LOG_D(MAC, "Allocating shared L1/L2 interface structure for instance %d @ %p\n", Mod_id, nr_if_inst[Mod_id]);

    nr_if_inst[Mod_id]->CC_mask=0;
    nr_if_inst[Mod_id]->NR_UL_indication = NR_UL_indication;
    if (NFAPI_MODE == NFAPI_MONOLITHIC)
      nr_if_inst[Mod_id]->NR_slot_indication = run_scheduler_monolithic;
    else if (NFAPI_MODE == NFAPI_MODE_PNF)
      nr_if_inst[Mod_id]->NR_slot_indication = pnf_send_slot_ind;
    else // NFAPI_MODE_VNF
      NULL;
    AssertFatal(pthread_mutex_init(&nr_if_inst[Mod_id]->if_mutex,NULL)==0,
                "allocation of nr_if_inst[%d]->if_mutex fails\n",Mod_id);
  }

  init_sched_response();

  return nr_if_inst[Mod_id];
}
