Sat Oct 23 18:07:03 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/SD2/A3/dat80.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2

$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr
$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       23 OCT 2021
Days until program expires : 176
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      800
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      700
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m80.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:   7.94996911480464        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9989E+02  1.9604E+02  3.6233E+02 -2.0364E+01  1.9787E+02  8.0102E+01 -1.1207E+02 -3.7929E+02 -5.9983E+01 -1.9033E+02
            -5.1632E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1919.96198919321        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1222E+00  9.6739E-01  7.8484E-01  1.1698E+00  8.3714E-01  7.7665E-01  9.8580E-01  1.0740E+00  8.8143E-01  9.5911E-01
             5.2946E+00
 PARAMETER:  2.1532E-01  6.6846E-02 -1.4228E-01  2.5680E-01 -7.7768E-02 -1.5276E-01  8.5700E-02  1.7141E-01 -2.6211E-02  5.8245E-02
             1.7667E+00
 GRADIENT:   1.5935E+02  2.3063E+01 -1.7493E+01  6.1903E+01 -1.1252E+01 -2.6512E+01  1.4528E+01  1.3204E+01  1.9238E+01  3.1072E+01
             5.9441E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2025.53649288220        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0516E+00  5.5307E-01  3.2428E-01  1.3068E+00  4.0351E-01  8.0959E-01  1.6551E+00  4.0389E-01  1.0980E+00  1.1338E-01
             4.3292E+00
 PARAMETER:  1.5029E-01 -4.9227E-01 -1.0261E+00  3.6761E-01 -8.0754E-01 -1.1123E-01  6.0385E-01 -8.0662E-01  1.9352E-01 -2.0770E+00
             1.5654E+00
 GRADIENT:   1.8426E+01  4.2965E+01 -8.3478E+01  1.7800E+02  1.0701E+02 -1.8005E+01  3.5920E+01  3.2040E+00  2.0822E+01  4.9269E-01
             5.0820E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2142.12464052473        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.4009E-01  3.0666E-01  1.8422E-01  1.2629E+00  2.5248E-01  8.8440E-01  1.7711E+00  5.1537E-01  1.3415E+00  8.7393E-02
             2.5475E+00
 PARAMETER:  3.8223E-02 -1.0820E+00 -1.5916E+00  3.3344E-01 -1.2764E+00 -2.2841E-02  6.7159E-01 -5.6287E-01  3.9382E-01 -2.3373E+00
             1.0351E+00
 GRADIENT:  -1.7437E+02  2.3295E+01 -1.4800E+02  2.3576E+02  2.5461E+02 -1.3615E+01  1.3841E+01 -2.1944E+01 -1.9643E+00 -7.5709E+00
            -2.7017E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2182.68033073179        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  9.7178E-01  4.0529E-01  3.1994E-01  1.1938E+00  3.6202E-01  8.9857E-01  1.8314E+00  5.4002E-01  9.8736E-01  1.4032E-01
             2.6541E+00
 PARAMETER:  7.1372E-02 -8.0315E-01 -1.0396E+00  2.7713E-01 -9.1604E-01 -6.9465E-03  7.0511E-01 -5.1614E-01  8.7282E-02 -1.8638E+00
             1.0761E+00
 GRADIENT:  -1.3754E+02 -4.0404E+00 -3.9769E+01  4.6701E+01  1.1139E+02  3.8049E+00  1.7003E+01 -5.6684E+00 -6.0862E+00 -2.6120E+00
             2.7001E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2210.95567069738        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.0149E+00  3.2738E-01  2.4164E-01  1.1534E+00  2.7691E-01  8.7875E-01  1.5011E+00  7.2131E-01  1.0868E+00  5.8837E-01
             2.4481E+00
 PARAMETER:  1.1478E-01 -1.0166E+00 -1.3203E+00  2.4270E-01 -1.1840E+00 -2.9254E-02  5.0618E-01 -2.2669E-01  1.8323E-01 -4.3039E-01
             9.9529E-01
 GRADIENT:  -2.0625E+01  3.3369E+01  5.5584E+00  5.2267E+01 -1.4321E+01 -4.8114E-01  1.3102E+01  1.1175E+01 -7.2614E+00  6.7271E+00
            -2.7930E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2220.48354152170        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0225E+00  2.6999E-01  1.9139E-01  1.0533E+00  2.3499E-01  8.7806E-01  1.3014E+00  2.7350E-01  1.2058E+00  6.5271E-01
             2.4861E+00
 PARAMETER:  1.2222E-01 -1.2094E+00 -1.5535E+00  1.5193E-01 -1.3482E+00 -3.0036E-02  3.6348E-01 -1.1964E+00  2.8718E-01 -3.2662E-01
             1.0107E+00
 GRADIENT:   8.2416E-01  3.0689E+00 -7.6516E+00 -9.4128E-01  8.9271E+00 -4.9438E-01  2.5502E+00  6.6515E-01  2.6508E+00  2.2840E+00
             2.8824E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2220.75032549532        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      902
 NPARAMETR:  1.0227E+00  2.6772E-01  1.9150E-01  1.0553E+00  2.3389E-01  8.7875E-01  1.2502E+00  1.1990E-01  1.1989E+00  6.6182E-01
             2.4888E+00
 PARAMETER:  1.2245E-01 -1.2178E+00 -1.5529E+00  1.5384E-01 -1.3529E+00 -2.9259E-02  3.2332E-01 -2.0211E+00  2.8137E-01 -3.1276E-01
             1.0118E+00
 GRADIENT:   1.1766E+00  7.2501E-01 -1.2256E+00  3.6052E+00 -5.4143E-01 -1.8422E-01 -3.2076E-01  6.1251E-02  1.1575E+00  8.9228E-01
            -1.2192E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2220.78556845380        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1077
 NPARAMETR:  1.0223E+00  2.6792E-01  1.9257E-01  1.0538E+00  2.3467E-01  8.7910E-01  1.2595E+00  2.8501E-02  1.1903E+00  6.5921E-01
             2.4925E+00
 PARAMETER:  1.2207E-01 -1.2171E+00 -1.5473E+00  1.5243E-01 -1.3496E+00 -2.8857E-02  3.3073E-01 -3.4578E+00  2.7422E-01 -3.1671E-01
             1.0133E+00
 GRADIENT:  -6.8750E-02 -3.4300E-02  1.3056E-01 -3.7453E-01  2.5232E-02  2.0631E-02  1.2746E-02  2.5221E-03 -4.2884E-03 -1.7695E-01
            -1.1545E-01

0ITERATION NO.:   44    OBJECTIVE VALUE:  -2220.78676460638        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1204
 NPARAMETR:  1.0223E+00  2.6792E-01  1.9252E-01  1.0540E+00  2.3463E-01  8.7906E-01  1.2584E+00  1.0000E-02  1.1905E+00  6.6006E-01
             2.4927E+00
 PARAMETER:  1.2210E-01 -1.2171E+00 -1.5475E+00  1.5264E-01 -1.3497E+00 -2.8903E-02  3.2986E-01 -4.6489E+00  2.7433E-01 -3.1543E-01
             1.0134E+00
 GRADIENT:   1.2293E-02  4.4482E-03 -1.7708E-02  5.1648E-02  5.5775E-03  7.2787E-04 -2.1519E-03  0.0000E+00 -1.4743E-03  7.8395E-03
            -9.6370E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1204
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.4501E-04  1.7575E-02  1.1838E-06 -6.7922E-03  9.6056E-03
 SE:             2.9313E-02  1.8801E-02  2.2765E-04  2.7732E-02  2.4168E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7700E-01  3.4987E-01  9.9585E-01  8.0652E-01  6.9104E-01

 ETASHRINKSD(%)  1.7972E+00  3.7015E+01  9.9237E+01  7.0944E+00  1.9033E+01
 ETASHRINKVR(%)  3.5620E+00  6.0329E+01  9.9994E+01  1.3685E+01  3.4444E+01
 EBVSHRINKSD(%)  1.9521E+00  3.7016E+01  9.9226E+01  5.6672E+00  1.9237E+01
 EBVSHRINKVR(%)  3.8662E+00  6.0330E+01  9.9994E+01  1.1013E+01  3.4774E+01
 RELATIVEINF(%)  9.6089E+01  1.1279E+01  5.7454E-04  4.4089E+01  4.5516E+00
 EPSSHRINKSD(%)  2.2299E+01
 EPSSHRINKVR(%)  3.9625E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          700
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1286.5139464865417     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2220.7867646063792     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -934.27281811983744     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2220.787       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.68E-01  1.93E-01  1.05E+00  2.35E-01  8.79E-01  1.26E+00  1.00E-02  1.19E+00  6.60E-01  2.49E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       73.473
Stop Time:
Sat Oct 23 18:07:17 CDT 2021
