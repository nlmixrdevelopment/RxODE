#pragma once
#ifndef __RxODE_control_H__
#define __RxODE_control_H__
#define Rxc_scale 0
#define Rxc_method 1
#define Rxc_transitAbs 2
#define Rxc_atol 3
#define Rxc_rtol 4
#define Rxc_maxsteps 5
#define Rxc_hmin 6
#define Rxc_hmax 7
#define Rxc_hini 8
#define Rxc_maxordn 9
#define Rxc_maxords 10
#define Rxc_covsInterpolation 11
#define Rxc_addCov 12
#define Rxc_matrix 13
#define Rxc_sigma 14
#define Rxc_sigmaDf 15
#define Rxc_nCoresRV 16
#define Rxc_sigmaIsChol 17
#define Rxc_sigmaSeparation 18
#define Rxc_sigmaXform 19
#define Rxc_nDisplayProgress 20
#define Rxc_amountUnits 21
#define Rxc_timeUnits 22
#define Rxc_theta 23
#define Rxc_eta 24
#define Rxc_addDosing 25
#define Rxc_stateTrim 26
#define Rxc_updateObject 27
#define Rxc_omega 28
#define Rxc_omegaDf 29
#define Rxc_omegaIsChol 30
#define Rxc_omegaSeparation 31
#define Rxc_omegaXform 32
#define Rxc_nSub 33
#define Rxc_thetaMat 34
#define Rxc_thetaDf 35
#define Rxc_thetaIsChol 36
#define Rxc_nStud 37
#define Rxc_dfSub 38
#define Rxc_dfObs 39
#define Rxc_seed 40
#define Rxc_nsim 41
#define Rxc_minSS 42
#define Rxc_maxSS 43
#define Rxc_strictSS 44
#define Rxc_infSSstep 45
#define Rxc_istateReset 46
#define Rxc_subsetNonmem 47
#define Rxc_hmaxSd 48
#define Rxc_maxAtolRtolFactor 49
#define Rxc_from 50
#define Rxc_to 51
#define Rxc_by 52
#define Rxc_length_out 53
#define Rxc_iCov 54
#define Rxc_keep 55
#define Rxc_keepF 56
#define Rxc_keepI 57
#define Rxc_drop 58
#define Rxc_warnDrop 59
#define Rxc_omegaLower 60
#define Rxc_omegaUpper 61
#define Rxc_sigmaLower 62
#define Rxc_sigmaUpper 63
#define Rxc_thetaLower 64
#define Rxc_thetaUpper 65
#define Rxc_indLinPhiM 66
#define Rxc_indLinPhiTol 67
#define Rxc_indLinMatExpType 68
#define Rxc_indLinMatExpOrder 69
#define Rxc_idFactor 70
#define Rxc_mxhnil 71
#define Rxc_hmxi 72
#define Rxc_warnIdSort 73
#define Rxc_ssAtol 74
#define Rxc_ssRtol 75
#define Rxc_safeZero 76
#define Rxc_cacheEvent 77
#define Rxc_sumType 78
#define Rxc_prodType 79
#define Rxc_sensType 80
#define Rxc_linDiff 81
#define Rxc_linDiffCentral 82
#define Rxc_resample 83
#define Rxc_resampleID 84
#define RxMv_params 0
#define RxMv_lhs 1
#define RxMv_state 2
#define RxMv_trans 3
#define RxMv_model 4
#define RxMv_ini 5
#define RxMv_podo 6
#define RxMv_dfdy 7
#define RxMv_sens 8
#define RxMv_state_ignore 9
#define RxMv_version 10
#define RxMv_normal_state 11
#define RxMv_needSort 12
#define RxMv_nMtime 13
#define RxMv_extraCmt 14
#define RxMv_stateExtra 15
#define RxMv_dvid 16
#define RxMv_indLin 17
#define RxMv_flags 18
#define RxMv_slhs 19
#define RxMv_timeId 20
#define RxMv_md5 21
#define RxMvFlag_ncmt 0
#define RxMvFlag_ka 1
#define RxMvFlag_linB 2
#define RxMvFlag_maxeta 3
#define RxMvFlag_maxtheta 4
#define RxMvFlag_hasCmt 5
#define RxMvFlag_linCmt 6
#define RxMvFlag_linCmtFlg 7

#define RxMvTrans_lib_name 0
#define RxMvTrans_jac 1
#define RxMvTrans_prefix 2
#define RxMvTrans_dydt 3
#define RxMvTrans_calc_jac 4
#define RxMvTrans_calc_lhs 5
#define RxMvTrans_model_vars 6
#define RxMvTrans_theta 7
#define RxMvTrans_inis 8
#define RxMvTrans_dydt_lsoda 9
#define RxMvTrans_calc_jac_lsoda 10
#define RxMvTrans_ode_solver_solvedata 11
#define RxMvTrans_ode_solver_get_solvedata 12
#define RxMvTrans_dydt_liblsoda 13
#define RxMvTrans_F 14
#define RxMvTrans_Lag 15
#define RxMvTrans_Rate 16
#define RxMvTrans_Dur 17
#define RxMvTrans_mtime 18
#define RxMvTrans_assignFuns 19
#define RxMvTrans_ME 20
#define RxMvTrans_IndF 21
#define RxTrans_ndose 0
#define RxTrans_nobs 1
#define RxTrans_nid 2
#define RxTrans_cov1 3
#define RxTrans_covParPos 4
#define RxTrans_covParPosTV 5
#define RxTrans_sub0 6
#define RxTrans_baseSize 7
#define RxTrans_nTv 8
#define RxTrans_lst 9
#define RxTrans_nme 10
#define RxTrans_covParPos0 11
#define RxTrans_covUnits 12
#define RxTrans_pars 13
#define RxTrans_allBolus 14
#define RxTrans_allInf 15
#define RxTrans_mxCmt 16
#define RxTrans_lib_name 17
#define RxTrans_addCmt 18
#define RxTrans_cmtInfo 19
#define RxTrans_idLvl 20
#define RxTrans_allTimeVar 21
#define RxTrans_keepDosingOnly 22
#define RxTrans_censAdd 23
#define RxTrans_limitAdd 24
#define RxTrans_levelInfo 25
#define RxTrans_idInfo 26
#define RxTransNames CharacterVector _en(27);_en[0]="ndose";_en[1]="nobs";_en[2]="nid";_en[3]="cov1";_en[4]="covParPos";_en[5]="covParPosTV";_en[6]="sub0";_en[7]="baseSize";_en[8]="nTv";_en[9]="lst";_en[10]="nme";_en[11]="covParPos0";_en[12]="covUnits";_en[13]="pars";_en[14]="allBolus";_en[15]="allInf";_en[16]="mxCmt";_en[17]="lib_name";_en[18]="addCmt";_en[19]="cmtInfo";_en[20]="idLvl";_en[21]="allTimeVar";_en[22]="keepDosingOnly";_en[23]="censAdd";_en[24]="limitAdd";_en[25]="levelInfo";_en[26]="idInfo";e.names() = _en;

#endif // __RxODE_control_H__
