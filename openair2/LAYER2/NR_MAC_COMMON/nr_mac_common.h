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

/*! \file mac.h
* \brief common MAC function prototypes
* \author Navid Nikaein and Raymond Knopp, WIE-TAI CHEN
* \date Dec. 2019
* \version 0.1
* \company Eurecom
* \email raymond.knopp@eurecom.fr
*/

#ifndef __LAYER2_NR_MAC_COMMON_H__
#define __LAYER2_NR_MAC_COMMON_H__

#include "NR_MIB.h"
#include "NR_CellGroupConfig.h"
#include "NR_UE-NR-Capability.h"
#include "nr_mac.h"
#include "common/utils/nr/nr_common.h"

#define NB_SRS_PERIOD         (18)
static const uint16_t srs_period[NB_SRS_PERIOD] = { 0, 1, 2, 4, 5, 8, 10, 16, 20, 32, 40, 64, 80, 160, 320, 640, 1280, 2560};

typedef enum {
  pusch_dmrs_pos0 = 0,
  pusch_dmrs_pos1 = 1,
  pusch_dmrs_pos2 = 2,
  pusch_dmrs_pos3 = 3,
} pusch_dmrs_AdditionalPosition_t;

typedef enum {
  pusch_len1 = 1,
  pusch_len2 = 2
} pusch_maxLength_t;

typedef struct {
  uint32_t format;
  uint32_t start_symbol;
  uint32_t N_t_slot;
  uint32_t N_dur;
  uint32_t N_RA_slot;
  uint32_t N_RA_sfn;
  uint32_t max_association_period;
  int x;
  int y;
  int y2;
  uint64_t s_map;
} nr_prach_info_t;

uint32_t get_Y(const NR_SearchSpace_t *ss, int slot, rnti_t rnti);

uint8_t get_BG(uint32_t A, uint16_t R);
uint32_t get_short_bsr_value(int idx);
uint32_t get_long_bsr_value(int idx);
int16_t fill_dmrs_mask(const NR_PDSCH_Config_t *pdsch_Config,
                       int dci_format,
                       int dmrs_TypeA_Position,
                       int NrOfSymbols,
                       int startSymbol,
                       mappingType_t mappingtype,
                       int length);
int get_slots_per_frame_from_scs(int scs);
bool is_ul_slot(const slot_t slot, const frame_structure_t *fs);
bool is_dl_slot(const slot_t slot, const frame_structure_t *fs);
bool is_mixed_slot(const slot_t slot, const frame_structure_t *fs);
int get_tdd_period_idx(NR_TDD_UL_DL_ConfigCommon_t *tdd);
void config_frame_structure(int mu,
                            NR_TDD_UL_DL_ConfigCommon_t *tdd_UL_DL_ConfigurationCommon,
                            uint8_t tdd_period,
                            uint8_t frame_type,
                            frame_structure_t *fs);

uint8_t compute_srs_resource_indicator(long *maxMIMO_Layers,
                                       NR_PUSCH_Config_t *pusch_Config,
                                       NR_SRS_Config_t *srs_config,
                                       nr_srs_feedback_t *srs_feedback,
                                       uint32_t *val);

uint8_t compute_precoding_information(NR_PUSCH_Config_t *pusch_Config,
                                      NR_SRS_Config_t *srs_config,
                                      dci_field_t srs_resource_indicator,
                                      nr_srs_feedback_t *srs_feedback,
                                      const uint8_t *nrOfLayers,
                                      int *tpmi,
                                      uint32_t *val);

NR_PDSCH_TimeDomainResourceAllocationList_t *get_dl_tdalist(const NR_UE_DL_BWP_t *DL_BWP,
                                                            int controlResourceSetId,
                                                            int ss_type,
                                                            nr_rnti_type_t rnti_type);

NR_PUSCH_TimeDomainResourceAllocationList_t *get_ul_tdalist(const NR_UE_UL_BWP_t *UL_BWP,
                                                            int controlResourceSetId,
                                                            int ss_type,
                                                            nr_rnti_type_t rnti_type);

NR_tda_info_t get_ul_tda_info(const NR_UE_UL_BWP_t *ul_bwp,
                              int controlResourceSetId,
                              int ss_type,
                              nr_rnti_type_t rnti_type,
                              int tda_index);

NR_tda_info_t get_dl_tda_info(const NR_UE_DL_BWP_t *dl_BWP,
                              int ss_type,
                              int tda_index,
                              int dmrs_typeA_pos,
                              int mux_pattern,
                              nr_rnti_type_t rnti_type,
                              int coresetid,
                              bool sib1);

uint8_t getRBGSize(uint16_t bwp_size, long rbg_size_config);

uint16_t nr_dci_size(const NR_UE_DL_BWP_t *DL_BWP,
                     const NR_UE_UL_BWP_t *UL_BWP,
                     const NR_UE_ServingCell_Info_t *sc_info,
                     long pdsch_HARQ_ACK_Codebook,
                     dci_pdu_rel15_t *dci_pdu,
                     nr_dci_format_t format,
                     nr_rnti_type_t rnti_type,
                     NR_ControlResourceSet_t *coreset,
                     int bwp_id,
                     int ss_type,
                     uint16_t cset0_bwp_size,
                     uint16_t alt_size);

uint16_t get_rb_bwp_dci(nr_dci_format_t format,
                        int ss_type,
                        uint16_t cset0_bwp_size,
                        uint16_t ul_bwp_size,
                        uint16_t dl_bwp_size,
                        uint16_t initial_ul_bwp_size,
                        uint16_t initial_dl_bwp_size);

void find_aggregation_candidates(uint8_t *aggregation_level,
                                 uint8_t *nr_of_candidates,
                                 const NR_SearchSpace_t *ss,
                                 int maxL);

uint16_t get_nr_prach_format_from_index(uint8_t index, uint32_t pointa, uint8_t unpaired);

bool get_nr_prach_sched_from_info(nr_prach_info_t info,
                                  int config_index,
                                  int frame,
                                  int slot,
                                  int mu,
                                  frequency_range_t freq_range,
                                  uint16_t *RA_sfn_index,
                                  uint8_t unpaired);

nr_prach_info_t get_nr_prach_occasion_info_from_index(uint8_t index, frequency_range_t freq_range, uint8_t unpaired);

uint8_t get_pusch_mcs_table(long *mcs_Table,
                            int is_tp,
                            int dci_format,
                            int rnti_type,
                            int target_ss,
                            bool config_grant);

uint8_t compute_nr_root_seq(NR_RACH_ConfigCommon_t *rach_config,
                            uint8_t nb_preambles,
                            uint8_t unpaired,
                            frequency_range_t);

int ul_ant_bits(NR_DMRS_UplinkConfig_t *NR_DMRS_UplinkConfig, long transformPrecoder);

uint8_t get_pdsch_mcs_table(long *mcs_Table, int dci_format, int rnti_type, int ss_type);

int get_format0(uint8_t index, uint8_t unpaired, frequency_range_t frequency_range);

uint16_t get_NCS(uint8_t index, uint16_t format, uint8_t restricted_set_config);
int compute_pucch_crc_size(int O_uci);
uint8_t get_l0_ul(uint8_t mapping_type, uint8_t dmrs_typeA_position);
int32_t get_l_prime(uint8_t duration_in_symbols, uint8_t mapping_type, pusch_dmrs_AdditionalPosition_t additional_pos, pusch_maxLength_t pusch_maxLength, uint8_t start_symbolt, uint8_t dmrs_typeA_position);

uint8_t get_L_ptrs(uint8_t mcs1, uint8_t mcs2, uint8_t mcs3, uint8_t I_mcs, uint8_t mcs_table);
uint8_t get_K_ptrs(uint32_t nrb0, uint32_t nrb1, uint32_t N_RB);

uint32_t nr_compute_tbs(uint16_t Qm,
                        uint16_t R,
			uint16_t nb_rb,
			uint16_t nb_symb_sch,
			uint16_t nb_dmrs_prb,
                        uint16_t nb_rb_oh,
                        uint8_t tb_scaling,
			uint8_t Nl);

/** \brief Computes Q based on I_MCS PDSCH and table_idx for downlink. Implements MCS Tables from 38.214. */
uint8_t nr_get_Qm_dl(uint8_t Imcs, uint8_t table_idx);
uint32_t nr_get_code_rate_dl(uint8_t Imcs, uint8_t table_idx);

/** \brief Computes Q based on I_MCS PDSCH and table_idx for uplink. Implements MCS Tables from 38.214. */
uint8_t nr_get_Qm_ul(uint8_t Imcs, uint8_t table_idx);
uint32_t nr_get_code_rate_ul(uint8_t Imcs, uint8_t table_idx);

uint16_t get_nr_srs_offset(NR_SRS_PeriodicityAndOffset_t periodicityAndOffset);
void get_monitoring_period_offset(const NR_SearchSpace_t *ss, int *period, int *offset);

uint32_t nr_compute_tbslbrm(uint16_t table,
			    uint16_t nb_rb,
		            uint8_t Nl);

void get_type0_PDCCH_CSS_config_parameters(NR_Type0_PDCCH_CSS_config_t *type0_PDCCH_CSS_config,
                                           frame_t frameP,
                                           const NR_MIB_t *mib,
                                           uint8_t num_slot_per_frame,
                                           uint8_t ssb_subcarrier_offset,
                                           uint16_t ssb_start_symbol,
                                           NR_SubcarrierSpacing_t scs_ssb,
                                           frequency_range_t frequency_range,
                                           int nr_band,
                                           uint32_t ssb_index,
                                           uint32_t ssb_period,
                                           uint32_t ssb_offset_point_a);

uint16_t get_ssb_start_symbol(const long band, NR_SubcarrierSpacing_t scs, int i_ssb);

NR_tda_info_t get_info_from_tda_tables(default_table_type_t table_type,
                                       int tda,
                                       int dmrs_TypeA_Position,
                                       int normal_CP);

NR_tda_info_t set_tda_info_from_list(NR_PDSCH_TimeDomainResourceAllocationList_t *tdalist, int tda_index);

default_table_type_t get_default_table_type(int mux_pattern);

void fill_coresetZero(NR_ControlResourceSet_t *coreset0, NR_Type0_PDCCH_CSS_config_t *type0_PDCCH_CSS_config);
void fill_searchSpaceZero(NR_SearchSpace_t *ss0,
                          int slots_per_frame,
                          NR_Type0_PDCCH_CSS_config_t *type0_PDCCH_CSS_config);

uint8_t get_pusch_nb_antenna_ports(NR_PUSCH_Config_t *pusch_Config,
                                   NR_SRS_Config_t *srs_config,
                                   dci_field_t srs_resource_indicator);

uint16_t compute_pucch_prb_size(uint8_t format,
                                uint8_t nr_prbs,
                                uint16_t O_csi,
                                uint16_t O_ack,
                                uint8_t O_sr,
                                NR_PUCCH_MaxCodeRate_t *maxCodeRate,
                                uint8_t Qm,
                                uint8_t n_symb,
                                uint8_t n_re_ctrl);

float get_max_code_rate(NR_PUCCH_MaxCodeRate_t *maxCodeRate);
int get_f3_dmrs_symbols(NR_PUCCH_Resource_t *pucchres, NR_PUCCH_Config_t *pucch_Config);

unsigned int get_delta_f_RA_long(const unsigned int format);
unsigned int get_N_RA_RB(const unsigned int delta_f_RA_PRACH, const unsigned int delta_f_PUSCH);

void find_period_offset_SR(const NR_SchedulingRequestResourceConfig_t *SchedulingReqRec, int *period, int *offset);

void csi_period_offset(NR_CSI_ReportConfig_t *csirep,
                       struct NR_CSI_ResourcePeriodicityAndOffset *periodicityAndOffset,
                       int *period, int *offset);

bool set_dl_ptrs_values(NR_PTRS_DownlinkConfig_t *ptrs_config,
                        uint16_t rbSize, uint8_t mcsIndex, uint8_t mcsTable,
                        uint8_t *K_ptrs, uint8_t *L_ptrs,uint8_t *portIndex,
                        uint8_t *nERatio,uint8_t *reOffset,
                        uint8_t NrOfSymbols);

bool set_ul_ptrs_values(NR_PTRS_UplinkConfig_t *ul_ptrs_config,
                        uint16_t rbSize,uint8_t mcsIndex, uint8_t mcsTable,
                        uint8_t *K_ptrs, uint8_t *L_ptrs,
                        uint8_t *reOffset, uint8_t *maxNumPorts, uint8_t *ulPower,
                        uint8_t NrOfSymbols);

/* \brief Set the transform precoding according to 6.1.3 of 3GPP TS 38.214 version 16.3.0 Release 16
@param    *current_UL_BWP  pointer to uplink bwp
@param    dci_format       dci format
@param    configuredGrant  indicates whether a configured grant was received or not
@returns                   transformPrecoding value */
long get_transformPrecoding(const NR_UE_UL_BWP_t *current_UL_BWP, nr_dci_format_t dci_format, uint8_t configuredGrant);

void compute_csi_bitlen(const NR_CSI_MeasConfig_t *csi_MeasConfig, nr_csi_report_t *csi_report_template);

uint16_t nr_get_csi_bitlen(nr_csi_report_t *csi_report);

uint32_t compute_PDU_length(uint32_t num_TLV, uint32_t total_length);

rnti_t nr_get_ra_rnti(uint8_t s_id, uint8_t t_id, uint8_t f_id, uint8_t ul_carrier_id);

rnti_t nr_get_MsgB_rnti(uint8_t s_id, uint8_t t_id, uint8_t f_id, uint8_t ul_carrier_id);

bool supported_bw_comparison(int bw_mhz, NR_SupportedBandwidth_t *supported_BW, long *support_90mhz);

int get_FeedbackDisabled(NR_DownlinkHARQ_FeedbackDisabled_r17_t *downlinkHARQ_FeedbackDisabled_r17, int harq_pid);

int get_nrofHARQ_ProcessesForPDSCH(const NR_UE_ServingCell_Info_t *sc_info);

int get_nrofHARQ_ProcessesForPUSCH(const NR_UE_ServingCell_Info_t *sc_info);

int nr_get_prach_or_ul_mu(const NR_MsgA_ConfigCommon_r16_t *msgacc, const NR_RACH_ConfigCommon_t *rach_ConfigCommon, const int ul_mu);

#endif
