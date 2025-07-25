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

/*
                                gnb_app.c
                             -------------------
  AUTHOR  : Laurent Winckel, Sebastien ROUX, Lionel GAUTHIER, Navid Nikaein, WEI-TAI CHEN
  COMPANY : EURECOM, NTUST
  EMAIL   : Lionel.Gauthier@eurecom.fr and Navid Nikaein, kroempa@gmail.com
*/

#include <string.h>
#include <stdio.h>
#include <nr_pdcp/nr_pdcp.h>
#include <softmodem-common.h>
#include <nr-softmodem.h>

#include "gnb_app.h"
#include "assertions.h"
#include "common/ran_context.h"

#include "common/utils/LOG/log.h"

#include "x2ap_eNB.h"
#include "intertask_interface.h"
#include "ngap_gNB.h"
#include "sctp_eNB_task.h"
#include "openair3/ocp-gtpu/gtp_itf.h"
#include "PHY/INIT/phy_init.h" 
#include "f1ap_cu_task.h"
#include "f1ap_du_task.h"
#include "nfapi/oai_integration/vendor_ext.h"
#include <openair2/LAYER2/nr_pdcp/nr_pdcp.h>
#include "openair2/LAYER2/nr_pdcp/nr_pdcp_oai_api.h"
#include "openair2/E1AP/e1ap.h"
#include "gnb_config.h"
#include "openair2/LAYER2/NR_MAC_gNB/mac_proto.h"

extern RAN_CONTEXT_t RC;

#define GNB_REGISTER_RETRY_DELAY 10

/*------------------------------------------------------------------------------*/


uint32_t gNB_app_register(uint32_t gnb_id_start, uint32_t gnb_id_end)//, const Enb_properties_array_t *enb_properties)
{
  uint32_t         gnb_id;
  MessageDef      *msg_p;
  uint32_t         register_gnb_pending = 0;

  for (gnb_id = gnb_id_start; (gnb_id < gnb_id_end) ; gnb_id++) {
    {
      if (IS_SA_MODE(get_softmodem_params())) {

        // note:  there is an implicit relationship between the data structure and the message name
        msg_p = itti_alloc_new_message (TASK_GNB_APP, 0, NGAP_REGISTER_GNB_REQ); //Message Temporarily reuse

        RCconfig_NR_NG(msg_p, gnb_id);

        itti_send_msg_to_task (TASK_NGAP, GNB_MODULE_ID_TO_INSTANCE(gnb_id), msg_p);
      }
    }

    register_gnb_pending++;
    }

  return register_gnb_pending;
}


/*------------------------------------------------------------------------------*/
uint32_t gNB_app_register_x2(uint32_t gnb_id_start, uint32_t gnb_id_end) {
  uint32_t         gnb_id;
  MessageDef      *msg_p;
  uint32_t         register_gnb_x2_pending = 0;

  for (gnb_id = gnb_id_start; (gnb_id < gnb_id_end) ; gnb_id++) {
    {
      msg_p = itti_alloc_new_message (TASK_GNB_APP, 0, X2AP_REGISTER_ENB_REQ);
      RCconfig_NR_X2(msg_p, gnb_id);
      itti_send_msg_to_task (TASK_X2AP, ENB_MODULE_ID_TO_INSTANCE(gnb_id), msg_p);
      register_gnb_x2_pending++;
    }
  }

  return register_gnb_x2_pending;
}

/*------------------------------------------------------------------------------*/

void *gNB_app_task(void *args_p)
{

  MessageDef                      *msg_p           = NULL;
  const char                      *msg_name        = NULL;
  instance_t                      instance;
  int                             result;
  /* for no gcc warnings */
  (void)instance;

  int cell_to_activate = 0;
  itti_mark_task_ready (TASK_GNB_APP);
  ngran_node_t node_type = get_node_type();

  if (RC.nb_nr_inst > 0) {
    if (node_type == ngran_gNB_CUCP ||
        node_type == ngran_gNB_CU ||
        node_type == ngran_eNB_CU ||
        node_type == ngran_ng_eNB_CU) {

      if (itti_create_task(TASK_CU_F1, F1AP_CU_task, NULL) < 0) {
        LOG_E(F1AP, "Create task for F1AP CU failed\n");
        AssertFatal(1==0,"exiting");
      }
    }

    if (node_type == ngran_gNB_CUCP) {
      if (itti_create_task(TASK_CUCP_E1, E1AP_CUCP_task, NULL) < 0)
        AssertFatal(false, "Create task for E1AP CP failed\n");
      E1_t e1type = CPtype;
      MessageDef *msg = RCconfig_NR_CU_E1(&e1type);
      AssertFatal(msg != NULL, "Send ITTI to task for E1AP CP failed\n");
      // this sends the E1AP_REGISTER_REQ to CU-CP so it sets up the socket
      // it does NOT use the E1AP part
      itti_send_msg_to_task(TASK_CUCP_E1, 0, msg);
    }

    if (node_type == ngran_gNB_CUUP) {
      AssertFatal(false, "To run CU-UP use executable nr-cuup\n");
    }

    if (NODE_IS_DU(node_type)) {
      if (itti_create_task(TASK_DU_F1, F1AP_DU_task, NULL) < 0) {
        LOG_E(F1AP, "Create task for F1AP DU failed\n");
        AssertFatal(1==0,"exiting");
      }
    }
    if (NODE_IS_DU(node_type) || NODE_IS_MONOLITHIC(node_type)) {
      // need to check SA?
      nr_mac_send_f1_setup_req();
    }
  }
  do {
    // Wait for a message
    itti_receive_msg (TASK_GNB_APP, &msg_p);

    msg_name = ITTI_MSG_NAME (msg_p);
    instance = ITTI_MSG_DESTINATION_INSTANCE (msg_p);

    switch (ITTI_MSG_ID(msg_p)) {
    case TERMINATE_MESSAGE:
      LOG_W(GNB_APP, " *** Exiting GNB_APP thread\n");
      itti_exit_task ();
      break;

    case MESSAGE_TEST:
      LOG_I(GNB_APP, "Received %s\n", ITTI_MSG_NAME(msg_p));
      break;



    case NGAP_REGISTER_GNB_CNF:
      LOG_I(GNB_APP, "[gNB %ld] Received %s: associated AMF %d\n", instance, msg_name,
            NGAP_REGISTER_GNB_CNF(msg_p).nb_amf);
      break;

    case F1AP_SETUP_RESP:
      AssertFatal(false, "Should not received this, logic bug\n");
      break;

    case F1AP_GNB_CU_CONFIGURATION_UPDATE:
      AssertFatal(NODE_IS_DU(node_type), "Should not have received F1AP_GNB_CU_CONFIGURATION_UPDATE in CU/gNB\n");
      LOG_I(GNB_APP,
            "Received %s: associated with %d cells to activate\n",
            ITTI_MSG_NAME(msg_p),
            F1AP_GNB_CU_CONFIGURATION_UPDATE(msg_p).num_cells_to_activate);
      cell_to_activate += F1AP_GNB_CU_CONFIGURATION_UPDATE(msg_p).num_cells_to_activate;
      gNB_app_handle_f1ap_gnb_cu_configuration_update(&F1AP_GNB_CU_CONFIGURATION_UPDATE(msg_p));

      /* Check if at least gNB is registered with one AMF */
      AssertFatal(cell_to_activate == 1,"No cells to activate or cells > 1 %d\n",cell_to_activate);

      break;

    case NGAP_DEREGISTERED_GNB_IND:
      LOG_W(GNB_APP, "[gNB %ld] Received %s: associated AMF %d\n", instance, msg_name,
            NGAP_DEREGISTERED_GNB_IND(msg_p).nb_amf);

      /* TODO handle recovering of registration */
      break;

    case TIMER_HAS_EXPIRED:
      LOG_I(GNB_APP, " Received %s: timer_id %ld\n", msg_name, TIMER_HAS_EXPIRED(msg_p).timer_id);
      break;

    default:
      LOG_E(GNB_APP, "Received unexpected message %s\n", msg_name);
      break;
    }

    result = itti_free (ITTI_MSG_ORIGIN_ID(msg_p), msg_p);
    AssertFatal (result == EXIT_SUCCESS, "Failed to free memory (%d)!\n", result);
  } while (1);


  return NULL;
}
