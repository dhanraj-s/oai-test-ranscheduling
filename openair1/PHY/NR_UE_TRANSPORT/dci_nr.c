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

/*! \file dci_nr.c
 * \brief Implements PDCCH physical channel TX/RX procedures (36.211) and DCI encoding/decoding (36.212/36.213). Current LTE
 * compliance V8.6 2009-03. \author R. Knopp, A. Mico Pereperez \date 2018 \version 0.1 \company Eurecom \email: knopp@eurecom.fr
 * \note
 * \warning
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "executables/softmodem-common.h"
#include "nr_transport_proto_ue.h"
#include "PHY/CODING/nrPolar_tools/nr_polar_dci_defs.h"
#include "PHY/phy_extern.h"
#include "PHY/CODING/coding_extern.h"
#include "PHY/sse_intrin.h"
#include "common/utils/nr/nr_common.h"
#include <openair1/PHY/TOOLS/phy_scope_interface.h>
#include "openair1/PHY/NR_REFSIG/nr_refsig_common.h"

#include "assertions.h"
#include "T.h"

//#define DEBUG_DCI_DECODING 1

//#define NR_PDCCH_DCI_DEBUG            // activates NR_PDCCH_DCI_DEBUG logs
#ifdef NR_PDCCH_DCI_DEBUG
#define LOG_DDD(a, ...) printf("<-NR_PDCCH_DCI_DEBUG (%s)-> " a, __func__, ##__VA_ARGS__ )
#define LOG_DSYMB(b)                                                               \
  LOG_DDD("RB[c_rb %d] \t RE[re %d] => rxF_ext[%d]=(%d,%d)\t rxF[%d]=(%d,%d)\n" b, \
          c_rb,                                                                    \
          i,                                                                       \
          j,                                                                       \
          rxF_ext[j].r,                                                            \
          rxF_ext[j].i,                                                            \
          i,                                                                       \
          rxF[i].r,                                                                \
          rxF[i].i)
#else
#define LOG_DDD(a...)
#define LOG_DSYMB(a...)
#endif
#define NR_NBR_CORESET_ACT_BWP 3 // The number of CoreSets per BWP is limited to 3 (including initial CORESET: ControlResourceId 0)
#define NR_NBR_SEARCHSPACE_ACT_BWP \
  10 // The number of SearSpaces per BWP is limited to 10 (including initial SEARCHSPACE: SearchSpaceId 0)

#ifdef LOG_I
#undef LOG_I
#define LOG_I(A, B...) printf(B)
#endif

static void nr_pdcch_demapping_deinterleaving(c16_t *llr,
                                              c16_t *e_rx,
                                              uint8_t coreset_time_dur,
                                              uint8_t start_symbol,
                                              uint32_t coreset_nbr_rb,
                                              uint8_t reg_bundle_size_L,
                                              uint8_t coreset_interleaver_size_R,
                                              uint8_t n_shift,
                                              uint8_t number_of_candidates,
                                              uint16_t *CCE,
                                              uint8_t *L)
{
  /*
   * This function will do demapping and deinterleaving from llr containing demodulated symbols
   * Demapping will regroup in REG and bundles
   * Deinterleaving will order the bundles
   *
   * In the following example we can see the process. The llr contains the demodulated IQs, but they are not ordered from
   REG 0,1,2,..
   * In e_rx (z) we will order the REG ids and group them into bundles.
   * Then we will put the bundles in the correct order as indicated in subclause 7.3.2.2
   *
   llr --------------------------> e_rx (z) ----> e_rx (z)
   |   ...
   |   ...
   |   REG 26
   symbol 2    |   ...
   |   ...
   |   REG 5
   |   REG 2

   |   ...
   |   ...
   |   REG 25
   symbol 1    |   ...
   |   ...
   |   REG 4
   |   REG 1

   |   ...
   |   ...                           ...              ...
   |   REG 24 (bundle 7)             ...              ...
   symbol 0    |   ...                           bundle 3         bundle 6
   |   ...                           bundle 2         bundle 1
   |   REG 3                         bundle 1         bundle 7
   |   REG 0  (bundle 0)             bundle 0         bundle 0

  */

  uint32_t coreset_C = 0;
  int coreset_interleaved = 0;
  const int N_regs = coreset_nbr_rb * coreset_time_dur;

  if (reg_bundle_size_L != 0) { // interleaving will be done only if reg_bundle_size_L != 0
    coreset_interleaved = 1;
    coreset_C = (uint32_t)(N_regs / (coreset_interleaver_size_R * reg_bundle_size_L));
  } else {
    reg_bundle_size_L = 6;
  }

  int B_rb = reg_bundle_size_L / coreset_time_dur; // nb of RBs occupied by each REG bundle
  int num_bundles_per_cce = 6 / reg_bundle_size_L;
  int n_cce = N_regs / 6;
  int max_bundles = n_cce * num_bundles_per_cce;
  int f_bundle_j_list[max_bundles];
  // for each bundle
  int c = 0, r = 0, f_bundle_j = 0;
  for (int nb = 0; nb < max_bundles; nb++) {
    if (coreset_interleaved == 0)
      f_bundle_j = nb;
    else {
      if (r == coreset_interleaver_size_R) {
        r = 0;
        c++;
      }
      f_bundle_j = ((r * coreset_C) + c + n_shift) % (N_regs / reg_bundle_size_L);
      r++;
    }
    f_bundle_j_list[nb] = f_bundle_j;
  }

  // Get cce_list indices by bundle index in ascending order
  int f_bundle_j_list_ord[number_of_candidates][max_bundles];
  for (int c_id = 0; c_id < number_of_candidates; c_id++) {
    int start_bund_cand = CCE[c_id] * num_bundles_per_cce;
    int max_bund_per_cand = L[c_id] * num_bundles_per_cce;
    int f_bundle_j_list_id = 0;
    for (int nb = 0; nb < max_bundles; nb++) {
      for (int bund_cand = start_bund_cand; bund_cand < start_bund_cand + max_bund_per_cand; bund_cand++) {
        if (f_bundle_j_list[bund_cand] == nb) {
          f_bundle_j_list_ord[c_id][f_bundle_j_list_id] = nb;
          f_bundle_j_list_id++;
        }
      }
    }
  }

  int rb_count = 0;
  const int data_sc = 9; // 9 sub-carriers with data per PRB
  for (int c_id = 0; c_id < number_of_candidates; c_id++) {
    for (int symbol_idx = start_symbol; symbol_idx < start_symbol + coreset_time_dur; symbol_idx++) {
      for (int cce_count = 0; cce_count < L[c_id]; cce_count++) {
        for (int k = 0; k < NR_NB_REG_PER_CCE / reg_bundle_size_L; k++) { // loop over REG bundles
          int f = f_bundle_j_list_ord[c_id][k + NR_NB_REG_PER_CCE * cce_count / reg_bundle_size_L];
          c16_t *in = llr + (f * B_rb + symbol_idx * coreset_nbr_rb) * data_sc;
          // loop over the RBs of the bundle
          memcpy(e_rx + data_sc * rb_count, in, B_rb * data_sc * sizeof(*e_rx));
          rb_count += B_rb;
        }
      }
    }
  }
}

static void nr_pdcch_llr(NR_DL_FRAME_PARMS *frame_parms,
                         int32_t rx_size,
                         c16_t rxdataF_comp[][rx_size],
                         c16_t *pdcch_llr,
                         uint8_t symbol,
                         uint32_t coreset_nbr_rb)
{
  c16_t *rxF = &rxdataF_comp[0][(symbol * coreset_nbr_rb * 12)];
  c16_t *pdcch_llrp = &pdcch_llr[symbol * coreset_nbr_rb * 9];

  if (!pdcch_llrp) {
    LOG_E(NR_PHY_DCI, "pdcch_qpsk_llr: llr is null, symbol %d\n", symbol);
    return;
  }

  LOG_DDD("llr logs: pdcch qpsk llr for symbol %d (pos %d), llr offset %ld\n",symbol,(symbol*frame_parms->N_RB_DL*12),pdcch_llrp-pdcch_llr);

  for (int i = 0; i < coreset_nbr_rb * 9; i++) {
    // We clip the signal
    c16_t res;
    res.r = min(rxF->r, 31);
    res.r = max(-32, res.r);
    res.i = min(rxF->i, 31);
    res.i = max(-32, res.i);
    *pdcch_llrp = res;
    LOG_DDD("llr logs: rb=%d i=%d *rxF:%d => *pdcch_llrp:%d\n",i/18,i,*rxF,*pdcch_llrp);
    rxF++;
    pdcch_llrp++;
  }
}

//compute average channel_level on each (TX,RX) antenna pair
void nr_pdcch_channel_level(int32_t rx_size,
                            c16_t dl_ch_estimates_ext[][rx_size],
                            NR_DL_FRAME_PARMS *frame_parms,
                            int32_t *avg,
                            int symbol,
                            int nb_rb)
{
  for (int aarx = 0; aarx < frame_parms->nb_antennas_rx; aarx++) {

    simde__m128i *dl_ch128 = (simde__m128i *)&dl_ch_estimates_ext[aarx][symbol * nb_rb * 12];

    //compute average level
    avg[aarx] = simde_mm_average(dl_ch128, nb_rb * 12, 0, nb_rb * 12);
    //LOG_DDD("Channel level : %d\n", avg[aarx]);
  }
}

// This function will extract the mapped DM-RS PDCCH REs as per 38.211 Section 7.4.1.3.2 (Mapping to physical resources)
static void nr_pdcch_extract_rbs_single(uint32_t rxdataF_sz,
                                        c16_t rxdataF[][rxdataF_sz],
                                        int32_t est_size,
                                        c16_t dl_ch_estimates[][est_size],
                                        int32_t rx_size,
                                        c16_t rxdataF_ext[][rx_size],
                                        c16_t dl_ch_estimates_ext[][rx_size],
                                        int symbol,
                                        NR_DL_FRAME_PARMS *frame_parms,
                                        uint8_t *coreset_freq_dom,
                                        uint32_t coreset_nbr_rb,
                                        uint32_t n_BWP_start)
{
  /*
   * This function is demapping DM-RS PDCCH RE
   * Implementing 38.211 Section 7.4.1.3.2 Mapping to physical resources
   * PDCCH DM-RS signals are mapped on RE a_k_l where:
   * k = 12*n + 4*kprime + 1
   * n=0,1,..
   * kprime=0,1,2
   * According to this equations, DM-RS PDCCH are mapped on k where k%12==1 || k%12==5 || k%12==9
   *
   */

#define NBR_RE_PER_RB_WITH_DMRS 12
  // after removing the 3 DMRS RE, the RB contains 9 RE with PDCCH
#define NBR_RE_PER_RB_WITHOUT_DMRS 9
  for (int aarx = 0; aarx < frame_parms->nb_antennas_rx; aarx++) {
    const c16_t *dl_ch0 = dl_ch_estimates[aarx] + frame_parms->ofdm_symbol_size * symbol;
    c16_t *rxFbase = rxdataF[aarx] + frame_parms->ofdm_symbol_size * symbol;
    LOG_DDD("dl_ch0 = &dl_ch_estimates[aarx = (%d)][0]\n", aarx);

    const int offset = symbol * coreset_nbr_rb * NBR_RE_PER_RB_WITH_DMRS;
    c16_t *dl_ch0_ext = &dl_ch_estimates_ext[aarx][offset];
    LOG_DDD("dl_ch0_ext = &dl_ch_estimates_ext[aarx = (%d)][symbol * (frame_parms->N_RB_DL * 9) = (%d)]\n", aarx, offset);
    c16_t *rxF_ext = &rxdataF_ext[aarx][offset];
    LOG_DDD("rxF_ext = &rxdataF_ext[aarx = (%d)][symbol * (frame_parms->N_RB_DL * 9) = (%d)]\n", aarx, offset);

    /*
     * The following for loop handles treatment of PDCCH contained in table rxdataF (in frequency domain)
     * In NR the PDCCH IQ symbols are contained within RBs in the CORESET defined by higher layers which is located within the BWP
     * Lets consider that the first RB to be considered as part of the CORESET and part of the PDCCH is n_BWP_start
     * Several cases have to be handled differently as IQ symbols are situated in different parts of rxdataF:
     * 1. Number of RBs in the system bandwidth is even
     *    1.1 The RB is <  than the N_RB_DL/2 -> IQ symbols are in the second half of the rxdataF (from first_carrier_offset)
     *    1.2 The RB is >= than the N_RB_DL/2 -> IQ symbols are in the first half of the rxdataF (from element 0)
     * 2. Number of RBs in the system bandwidth is odd
     * (particular case when the RB with DC as it is treated differently: it is situated in symbol borders of rxdataF)
     *    2.1 The RB is <  than the N_RB_DL/2 -> IQ symbols are in the second half of the rxdataF (from first_carrier_offset)
     *    2.2 The RB is >  than the N_RB_DL/2 -> IQ symbols are in the first half of the rxdataF (from element 0 + 2nd half RB
     * containing DC) 2.3 The RB is == N_RB_DL/2          -> IQ symbols are in the upper border of the rxdataF for first 6 IQ
     * element and the lower border of the rxdataF for the last 6 IQ elements If the first RB containing PDCCH within the UE BWP
     * and within the CORESET is higher than half of the system bandwidth (N_RB_DL), then the IQ symbol is going to be found at
     * the position 0+c_rb-N_RB_DL/2 in rxdataF and we have to point the pointer at (1+c_rb-N_RB_DL/2) in rxdataF
     */

    int c_rb = 0;
    for (int rb = 0; rb < coreset_nbr_rb; rb++, c_rb++) {
      int c_rb_by6 = c_rb / 6;

      // skip zeros in frequency domain bitmap
      while ((coreset_freq_dom[c_rb_by6 >> 3] & (1 << (7 - (c_rb_by6 & 7)))) == 0) {
        c_rb += 6;
        c_rb_by6 = c_rb / 6;
      }

      // first we set initial conditions for pointer to rxdataF depending on the situation of the first RB within the CORESET
      // (c_rb = n_BWP_start)
      if ((frame_parms->N_RB_DL & 1) == 1 && (c_rb + n_BWP_start) == (frame_parms->N_RB_DL >> 1)) {
        // treatment of RB containing the DC
        // if odd number RBs in system bandwidth and first RB to be treated is higher than middle system bandwidth (around DC)
        // we have to treat the RB in two parts: first part from i=0 to 5, the data is at the end of rxdataF (pointing at the
        // end of the table)
        c16_t *rxF = rxFbase + frame_parms->first_carrier_offset + 12 * (c_rb + n_BWP_start);

        int i = 0, j = 0;
        for (; i < 6; i++) { // treating first part of the RB note that i=5 would correspond to DC. We treat it in NR
          if ((i != 1) && (i != 5)) {
            dl_ch0_ext[j] = dl_ch0[i];
            rxF_ext[j] = rxF[i];
            LOG_DSYMB("");
            j++;
          } else {
            LOG_DSYMB("\t\t <==> DM-RS PDCCH, this is a pilot symbol\n");
          }
        }

        // then we point at the begining of the symbol part of rxdataF do process second part of RB
        for (; i < 12; i++) {
          if ((i != 9)) {
            dl_ch0_ext[j] = dl_ch0[i];
            rxF_ext[j] = rxFbase[i - 6];
            LOG_DSYMB("");
            j++;
          } else {
            LOG_DSYMB("\t\t <==> DM-RS PDCCH, this is a pilot symbol\n");
          }
        }

      } else { // treatment of any RB that does not contain the DC
        c16_t *rxF=NULL;
        if ((frame_parms->N_RB_DL & 1) == 0) {
          if ((c_rb + n_BWP_start) < (frame_parms->N_RB_DL >> 1))
            // if RB to be treated is lower than middle system bandwidth then rxdataF pointed
            // at (offset + c_br + symbol * ofdm_symbol_size): even case
            rxF = rxFbase + (frame_parms->first_carrier_offset + 12 * c_rb) + n_BWP_start * 12;
          else
            // number of RBs is even  and c_rb is higher than half system bandwidth (we don't skip DC)
            // if these conditions are true the pointer has to be situated at the 1st part of the rxdataF
            // we point at the 1st part of the rxdataF in symbol
            rxF = rxFbase + 12 * (c_rb + n_BWP_start - (frame_parms->N_RB_DL >> 1));
        } else {
          if ((c_rb + n_BWP_start) < (frame_parms->N_RB_DL >> 1))
            // if RB to be treated is lower than middle system bandwidth then rxdataF pointed
            //  at (offset + c_br + symbol * ofdm_symbol_size): odd case
            rxF = rxFbase + frame_parms->first_carrier_offset + 12 * (c_rb + n_BWP_start);
          else if ((c_rb + n_BWP_start) > (frame_parms->N_RB_DL >> 1))
            // number of RBs is odd  and   c_rb is higher than half system bandwidth + 1
            // if these conditions are true the pointer has to be situated at the 1st part of
            // the rxdataF just after the first IQ symbols of the RB containing DC
            // we point at the 1st part of the rxdataF in symbol
            rxF = rxFbase + 12 * (c_rb + n_BWP_start - (frame_parms->N_RB_DL >> 1)) - 6;
        }
	AssertFatal(rxF, "bug");
        int j = 0;

        for (int i = 0; i < 12; i++) {
          if ((i != 1) && (i != 5) && (i != 9)) {
            rxF_ext[j] = rxF[i];
            dl_ch0_ext[j] = dl_ch0[i];
            LOG_DSYMB("");
            j++;
          } else {
            LOG_DSYMB("\t\t <==> DM-RS PDCCH, this is a pilot symbol\n");
          }
        }
      }
      dl_ch0_ext += NBR_RE_PER_RB_WITHOUT_DMRS;
      rxF_ext += NBR_RE_PER_RB_WITHOUT_DMRS;
      dl_ch0 += 12;
    }
  }
}

#define print_shorts(s,x) printf("%s %d,%d,%d,%d,%d,%d,%d,%d\n",s,(x)[0],(x)[1],(x)[2],(x)[3],(x)[4],(x)[5],(x)[6],(x)[7])

void nr_pdcch_channel_compensation(int32_t rx_size,
                                   c16_t rxdataF_ext[][rx_size],
                                   c16_t dl_ch_estimates_ext[][rx_size],
                                   c16_t rxdataF_comp[][rx_size],
                                   int32_t **rho,
                                   NR_DL_FRAME_PARMS *frame_parms,
                                   uint8_t symbol,
                                   uint8_t output_shift,
                                   uint32_t coreset_nbr_rb)
{
  for (int aarx = 0; aarx < frame_parms->nb_antennas_rx; aarx++) {
    simde__m128i *dl_ch128 = (simde__m128i *)&dl_ch_estimates_ext[aarx][symbol * coreset_nbr_rb * 12];
    simde__m128i *rxdataF128 = (simde__m128i *)&rxdataF_ext[aarx][symbol * coreset_nbr_rb * 12];
    simde__m128i *rxdataF_comp128 = (simde__m128i *)&rxdataF_comp[aarx][symbol * coreset_nbr_rb * 12];
    //printf("ch compensation dl_ch ext addr %p \n", &dl_ch_estimates_ext[(aatx<<1)+aarx][symbol*20*12]);
    //printf("rxdataf ext addr %p symbol %d\n", &rxdataF_ext[aarx][symbol*20*12], symbol);
    //printf("rxdataf_comp addr %p\n",&rxdataF_comp[(aatx<<1)+aarx][symbol*20*12]);
    
    // multiply by conjugated channel
    mult_cpx_conj_vector((c16_t *)dl_ch128, (c16_t *)rxdataF128, (c16_t *)rxdataF_comp128, 12 * coreset_nbr_rb, output_shift);

    for (int rb = 0; rb < 12 * coreset_nbr_rb; rb++) {
        LOG_DDD("rxdataF128[%d]=(%d,%d) X dlch[%d]=(%d,%d) rxdataF_comp128[%d]=(%d,%d)\n",
                rb,
                ((c16_t *)rxdataF128)[rb].r,
                ((c16_t *)rxdataF128)[rb].i,
                rb,
                ((c16_t *)dl_ch128)[rb].r,
                ((c16_t *)dl_ch128)[rb].i,
                rb,
                ((c16_t *)rxdataF_comp128)[rb].r,
                ((c16_t *)rxdataF_comp128)[rb].i);

    }
  }
}

static void nr_pdcch_detection_mrc(NR_DL_FRAME_PARMS *frame_parms, int32_t rx_size, c16_t rxdataF_comp[][rx_size], int symbol)
{
  if (frame_parms->nb_antennas_rx>1) {
    const int sz = frame_parms->N_RB_DL * 12;
    
    simde__m128i *rxdataF_comp128_0 = (simde__m128i *)&rxdataF_comp[0][symbol * sz];
    simde__m128i *rxdataF_comp128_1 = (simde__m128i *)&rxdataF_comp[1][symbol * sz];

    // MRC on each re of rb
    for (int i = 0; i < sz >> 2; i++) {
      rxdataF_comp128_0[i] = simde_mm_adds_epi16(simde_mm_srai_epi16(rxdataF_comp128_0[i], 1), simde_mm_srai_epi16(rxdataF_comp128_1[i], 1));
    }
  }
}

void nr_rx_pdcch(PHY_VARS_NR_UE *ue,
                 const UE_nr_rxtx_proc_t *proc,
                 int32_t pdcch_est_size,
                 c16_t pdcch_dl_ch_estimates[][pdcch_est_size],
                 c16_t *pdcch_e_rx,
                 fapi_nr_dl_config_dci_dl_pdu_rel15_t *rel15,
                 c16_t rxdataF[][ue->frame_parms.samples_per_slot_wCP])
{
  NR_DL_FRAME_PARMS *frame_parms = &ue->frame_parms;

  uint8_t log2_maxh, aarx;
  int32_t avgs;
  int32_t avgP[4];
  int n_rb,rb_offset;
  get_coreset_rballoc(rel15->coreset.frequency_domain_resource,&n_rb,&rb_offset);

  // Pointers to extracted PDCCH symbols in frequency-domain.
  int32_t rx_size = ((4 * frame_parms->N_RB_DL * 12 + 31) >> 5) << 5;
  __attribute__((aligned(32))) c16_t rxdataF_ext[frame_parms->nb_antennas_rx][rx_size];
  __attribute__((aligned(32))) c16_t rxdataF_comp[frame_parms->nb_antennas_rx][rx_size];
  __attribute__((aligned(32))) c16_t pdcch_dl_ch_estimates_ext[frame_parms->nb_antennas_rx][rx_size];

  memset(rxdataF_comp, 0, sizeof(rxdataF_comp));

  // Pointer to llrs, 4-bit resolution.
  int32_t llr_size = 4 * n_rb * 9;
  c16_t llr[llr_size];

  memset(llr, 0, sizeof(llr));

  LOG_D(NR_PHY_DCI,
        "pdcch coreset: freq %x, n_rb %d, rb_offset %d\n",
        rel15->coreset.frequency_domain_resource[0],
        n_rb,
        rb_offset);
  for (int s=rel15->coreset.StartSymbolIndex; s<(rel15->coreset.StartSymbolIndex+rel15->coreset.duration); s++) {
    LOG_D(NR_PHY_DCI, "in nr_pdcch_extract_rbs_single(rxdataF -> rxdataF_ext || dl_ch_estimates -> dl_ch_estimates_ext)\n");

    nr_pdcch_extract_rbs_single(ue->frame_parms.samples_per_slot_wCP,
                                rxdataF,
                                pdcch_est_size,
                                pdcch_dl_ch_estimates,
                                rx_size,
                                rxdataF_ext,
                                pdcch_dl_ch_estimates_ext,
                                s,
                                frame_parms,
                                rel15->coreset.frequency_domain_resource,
                                n_rb,
                                rel15->BWPStart);

    LOG_D(NR_PHY_DCI,
          "we enter nr_pdcch_channel_level(avgP=%d) => compute channel level based on ofdm symbol 0, "
          "pdcch_vars[eNB_id]->dl_ch_estimates_ext\n",
          *avgP);
    LOG_D(NR_PHY_DCI, "in nr_pdcch_channel_level(dl_ch_estimates_ext -> dl_ch_estimates_ext)\n");
    // compute channel level based on ofdm symbol 0
    nr_pdcch_channel_level(rx_size,
                           pdcch_dl_ch_estimates_ext,
                           frame_parms,
                           avgP,
                           s,
                           n_rb);
    avgs = 0;

    for (aarx = 0; aarx < frame_parms->nb_antennas_rx; aarx++)
      avgs = cmax(avgs, avgP[aarx]);

    log2_maxh = (log2_approx(avgs) / 2) + 5;  //+frame_parms->nb_antennas_rx;

#ifdef UE_DEBUG_TRACE
    LOG_D(NR_PHY_DCI, "slot %d: pdcch log2_maxh = %d (%d,%d)\n", proc->nr_slot_rx, log2_maxh, avgP[0], avgs);
#endif
#if T_TRACER
    T(T_UE_PHY_PDCCH_ENERGY, T_INT(0), T_INT(0), T_INT(proc->frame_rx % 1024), T_INT(proc->nr_slot_rx), T_INT(avgP[0]), T_INT(avgP[1]), T_INT(avgP[2]), T_INT(avgP[3]));
#endif
    LOG_D(NR_PHY_DCI, "we enter nr_pdcch_channel_compensation(log2_maxh=%d)\n", log2_maxh);
    LOG_D(NR_PHY_DCI, "in nr_pdcch_channel_compensation(rxdataF_ext x dl_ch_estimates_ext -> rxdataF_comp)\n");
    // compute LLRs for ofdm symbol 0 only
    nr_pdcch_channel_compensation(rx_size, rxdataF_ext,
                                  pdcch_dl_ch_estimates_ext,
                                  rxdataF_comp,
                                  NULL,
                                  frame_parms,
                                  s,
                                  log2_maxh,
                                  n_rb); // log2_maxh+I0_shift

    UEscopeCopy(ue, pdcchRxdataF_comp, rxdataF_comp, sizeof(struct complex16), frame_parms->nb_antennas_rx, rx_size, 0);

    if (frame_parms->nb_antennas_rx > 1) {
      LOG_D(NR_PHY_DCI, "we enter nr_pdcch_detection_mrc(frame_parms->nb_antennas_rx=%d)\n", frame_parms->nb_antennas_rx);
      nr_pdcch_detection_mrc(frame_parms, rx_size, rxdataF_comp,s);
    }

    LOG_D(NR_PHY_DCI, "we enter nr_pdcch_llr(for symbol %d), pdcch_vars[eNB_id]->rxdataF_comp ---> pdcch_vars[eNB_id]->llr \n", s);
    LOG_D(NR_PHY_DCI, "in nr_pdcch_llr(rxdataF_comp -> llr)\n");
    nr_pdcch_llr(frame_parms,
                 rx_size,
                 rxdataF_comp,
                 llr,
                 s,
                 n_rb);

    UEscopeCopy(ue, pdcchLlr, llr, sizeof(c16_t), 1, llr_size, 0);

#ifdef DEBUG_DCI_DECODING
    printf("demapping: slot %u, mi %d\n",slot,get_mi(frame_parms,slot));
#endif
  }

  nr_pdcch_demapping_deinterleaving(llr,
                                    pdcch_e_rx,
                                    rel15->coreset.duration,
                                    rel15->coreset.StartSymbolIndex,
                                    n_rb,
                                    rel15->coreset.RegBundleSize,
                                    rel15->coreset.InterleaverSize,
                                    rel15->coreset.ShiftIndex,
                                    rel15->number_of_candidates,
                                    rel15->CCE,
                                    rel15->L);
}

static void nr_pdcch_unscrambling(c16_t *e_rx,
                                  uint16_t scrambling_RNTI,
                                  uint32_t length,
                                  uint16_t pdcch_DMRS_scrambling_id,
                                  int16_t *z2)
{
  uint32_t rnti = (uint32_t) scrambling_RNTI;
  uint16_t n_id = pdcch_DMRS_scrambling_id;
  uint32_t *seq = gold_cache(((rnti << 16) + n_id) % (1U << 31), length / 32); // this is c_init in 38.211 v15.1.0 Section 7.3.2.3
  LOG_D(NR_PHY_DCI, "PDCCH Unscrambling: scrambling_RNTI %x\n", rnti);
  int16_t *ptr = &e_rx[0].r;
  for (int i = 0; i < length; i++) {
    if (seq[i / 32] & (1UL << (i % 32)))
      z2[i] = -ptr[i];
    else
      z2[i] = ptr[i];
  }
}

void nr_dci_decoding_procedure(PHY_VARS_NR_UE *ue,
                               const UE_nr_rxtx_proc_t *proc,
                               c16_t *pdcch_e_rx,
                               fapi_nr_dci_indication_t *dci_ind,
                               fapi_nr_dl_config_dci_dl_pdu_rel15_t *rel15)
{
  int e_rx_cand_idx = 0;
  *dci_ind = (fapi_nr_dci_indication_t){.SFN = proc->frame_rx, .slot = proc->nr_slot_rx};
  // if DCI for SIB we don't break after finding 1st DCI with that RNTI
  // there might be SIB1 and otherSIB in the same slot with the same length
  bool is_SI = rel15->rnti == SI_RNTI;

  for (int j = 0; j < rel15->number_of_candidates; j++) {
    int CCEind = rel15->CCE[j];
    int L = rel15->L[j];

    // Loop over possible DCI lengths
    
    for (int k = 0; k < rel15->num_dci_options; k++) {
      // skip this candidate if we've already found one with the
      // same rnti and size at a different aggregation level
      int dci_length = rel15->dci_length_options[k];
      int ind;
      for (ind = 0; ind < dci_ind->number_of_dcis; ind++) {
        if (!is_SI && rel15->rnti == dci_ind->dci_list[ind].rnti && dci_length == dci_ind->dci_list[ind].payloadSize) {
          break;
        }
      }
      if (ind < dci_ind->number_of_dcis)
        continue;

      uint64_t dci_estimation[2] = {0};
      LOG_D(NR_PHY_DCI,
            "(%i.%i) Trying DCI candidate %d of %d number of candidates, CCE %d (%d), L %d, length %d, format %d\n",
            proc->frame_rx,
            proc->nr_slot_rx,
            j,
            rel15->number_of_candidates,
            CCEind,
            e_rx_cand_idx,
            L,
            dci_length,
            rel15->dci_format_options[k]);

      int16_t tmp_e[16 * 108];
      nr_pdcch_unscrambling(&pdcch_e_rx[e_rx_cand_idx],
                            rel15->coreset.scrambling_rnti,
                            L * 108,
                            rel15->coreset.pdcch_dmrs_scrambling_id,
                            tmp_e);

      const uint32_t crc = polar_decoder_int16(tmp_e, dci_estimation, 1, NR_POLAR_DCI_MESSAGE_TYPE, dci_length, L);

      rnti_t n_rnti = rel15->rnti;
      if (crc == n_rnti) {
        LOG_D(NR_PHY_DCI,
              "(%i.%i) Received dci indication (rnti %x,dci format %d,n_CCE %d,payloadSize %d,payload %llx)\n",
              proc->frame_rx,
              proc->nr_slot_rx,
              n_rnti,
              rel15->dci_format_options[k],
              CCEind,
              dci_length,
              *(unsigned long long *)dci_estimation);
        AssertFatal(dci_ind->number_of_dcis < sizeofArray(dci_ind->dci_list), "Fix allocation\n");
        fapi_nr_dci_indication_pdu_t *dci = dci_ind->dci_list + dci_ind->number_of_dcis;
        *dci = (fapi_nr_dci_indication_pdu_t){
            .rnti = n_rnti,
            .n_CCE = CCEind,
            .N_CCE = L,
            .dci_format = rel15->dci_format_options[k],
            .ss_type = rel15->ss_type_options[k],
            .coreset_type = rel15->coreset.CoreSetType,
        };
        int n_rb, rb_offset;
        get_coreset_rballoc(rel15->coreset.frequency_domain_resource, &n_rb, &rb_offset);
        dci->cset_start = rel15->BWPStart + rb_offset;
        dci->payloadSize = dci_length;
        memcpy(dci->payloadBits, dci_estimation, (dci_length + 7) / 8);
        dci_ind->number_of_dcis++;
        break;    // If DCI is found, no need to check for remaining DCI lengths
      } else {
        LOG_D(NR_PHY_DCI,
              "(%i.%i) Decoded crc %x does not match rnti %x for DCI format %d\n",
              proc->frame_rx,
              proc->nr_slot_rx,
              crc,
              n_rnti,
              rel15->dci_format_options[k]);
      }
    }
    e_rx_cand_idx += 9 * L * 6; // e_rx index for next candidate (L CCEs, 6 REGs per CCE and 9 REs per REG )
  }
}
