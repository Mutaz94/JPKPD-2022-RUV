Sun Oct 24 04:23:38 CDT 2021
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
$DATA ../../../../data/SD4/D/dat71.csv ignore=@
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
Current Date:       24 OCT 2021
Days until program expires : 175
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
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
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
 RAW OUTPUT FILE (FILE): m71.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1580.58907403420        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7614E+02 -1.2641E+01  6.1149E+01 -6.1152E+01 -1.2188E+01  3.1070E+01 -2.9735E+01 -1.5910E+01 -2.9169E+01 -2.8736E+01
             2.5336E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1593.54979630149        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.5151E-01  9.6884E-01  8.0539E-01  1.0594E+00  9.1646E-01  1.0630E+00  1.2671E+00  1.0858E+00  1.1684E+00  1.2334E+00
             8.4117E-01
 PARAMETER:  5.0293E-02  6.8342E-02 -1.1643E-01  1.5774E-01  1.2766E-02  1.6106E-01  3.3674E-01  1.8230E-01  2.5565E-01  3.0974E-01
            -7.2958E-02
 GRADIENT:  -8.3110E+00 -3.2615E+01 -2.8411E+01 -1.3558E+01  1.8996E+00  1.3089E+01 -3.2908E+00  1.4787E+01  2.8178E+01  2.2920E+01
            -2.3150E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1598.26800882535        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  9.4900E-01  7.4910E-01  8.0926E-01  1.2498E+00  8.1426E-01  1.0740E+00  1.8410E+00  7.3891E-01  9.4669E-01  1.1858E+00
             8.4248E-01
 PARAMETER:  4.7655E-02 -1.8889E-01 -1.1163E-01  3.2302E-01 -1.0547E-01  1.7136E-01  7.1029E-01 -2.0259E-01  4.5218E-02  2.7044E-01
            -7.1400E-02
 GRADIENT:  -1.0080E+01  1.4996E+01 -2.9063E+01  6.1260E+01  6.1198E+00  1.7740E+01  5.5939E+00  6.8489E+00  1.1806E+01  1.8029E+01
            -2.2201E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1603.06712039146        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.5630E-01  8.4262E-01  7.2194E-01  1.1332E+00  7.8811E-01  1.0142E+00  1.6970E+00  6.1473E-01  9.0926E-01  9.7330E-01
             8.9436E-01
 PARAMETER:  5.5321E-02 -7.1236E-02 -2.2581E-01  2.2504E-01 -1.3812E-01  1.1408E-01  6.2889E-01 -3.8658E-01  4.8758E-03  7.2934E-02
            -1.1647E-02
 GRADIENT:   4.0597E-01 -3.5028E-01  1.5713E+00 -7.7668E+00 -6.0962E+00 -4.8981E+00  3.3268E+00  2.2865E+00  8.6460E-01  1.8632E+00
             1.2301E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1603.76601084819        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  9.5840E-01  1.0097E+00  5.9212E-01  1.0275E+00  7.7248E-01  1.0297E+00  1.4478E+00  2.5765E-01  9.5572E-01  9.1752E-01
             8.9406E-01
 PARAMETER:  5.7507E-02  1.0970E-01 -4.2405E-01  1.2710E-01 -1.5815E-01  1.2923E-01  4.7004E-01 -1.2561E+00  5.4710E-02  1.3920E-02
            -1.1981E-02
 GRADIENT:   9.7675E-01  6.0311E+00  5.5308E-01  4.4156E+00 -5.0750E+00  2.6454E-01  1.7896E-01  2.2216E-01 -4.3714E-01  8.9684E-01
             3.8500E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1603.82338290506        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  9.5785E-01  1.0015E+00  5.9386E-01  1.0267E+00  7.7505E-01  1.0288E+00  1.4546E+00  2.1169E-01  9.5769E-01  9.1842E-01
             8.9360E-01
 PARAMETER:  5.6933E-02  1.0153E-01 -4.2112E-01  1.2638E-01 -1.5483E-01  1.2835E-01  4.7476E-01 -1.4526E+00  5.6764E-02  1.4896E-02
            -1.2494E-02
 GRADIENT:  -5.2815E-02  1.1804E-01 -1.6619E-02 -2.4946E-01  2.9166E-02 -7.0368E-02 -2.8585E-02  4.6299E-02  1.0334E-01 -5.1018E-02
            -7.9621E-03

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1603.82951992686        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  9.5791E-01  9.9942E-01  5.8574E-01  1.0265E+00  7.6823E-01  1.0290E+00  1.4558E+00  1.5116E-01  9.5566E-01  9.1228E-01
             8.9382E-01
 PARAMETER:  5.6996E-02  9.9417E-02 -4.3487E-01  1.2613E-01 -1.6367E-01  1.2857E-01  4.7554E-01 -1.7894E+00  5.4652E-02  8.1887E-03
            -1.2253E-02
 GRADIENT:  -1.7485E-02 -3.2020E-02 -6.6853E-02  3.3228E-03  3.7852E-02  3.1431E-03  2.4485E-02  2.3238E-03 -4.8810E-03  2.0419E-02
             7.9536E-03

0ITERATION NO.:   34    OBJECTIVE VALUE:  -1603.83026452500        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:     1184
 NPARAMETR:  9.5888E-01  1.0006E+00  5.8498E-01  1.0255E+00  7.6818E-01  1.0300E+00  1.4556E+00  1.4722E-01  9.5616E-01  9.1181E-01
             8.9382E-01
 PARAMETER:  5.8013E-02  1.0065E-01 -4.3619E-01  1.2520E-01 -1.6373E-01  1.2952E-01  4.7540E-01 -1.8158E+00  5.5171E-02  7.6713E-03
            -1.2251E-02
 GRADIENT:   5.3521E-01  2.2497E-02  4.6978E-02 -4.0095E-02 -1.0692E-01  6.2759E-02  1.5135E-02  1.7573E-04  2.8383E-02  4.4680E-02
             8.3270E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1184
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.9248E-05 -6.0633E-04 -8.4332E-03 -1.2995E-03 -9.9381E-03
 SE:             2.9879E-02  2.3261E-02  3.0538E-03  2.5011E-02  2.2083E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9735E-01  9.7920E-01  5.7528E-03  9.5856E-01  6.5268E-01

 ETASHRINKSD(%)  1.0000E-10  2.2073E+01  8.9769E+01  1.6210E+01  2.6021E+01
 ETASHRINKVR(%)  1.0000E-10  3.9274E+01  9.8953E+01  2.9793E+01  4.5270E+01
 EBVSHRINKSD(%)  3.2353E-01  2.1005E+01  9.1580E+01  1.6956E+01  2.5433E+01
 EBVSHRINKVR(%)  6.4601E-01  3.7599E+01  9.9291E+01  3.1038E+01  4.4397E+01
 RELATIVEINF(%)  9.9275E+01  6.7711E+00  1.1041E-01  8.5295E+00  5.5388E+00
 EPSSHRINKSD(%)  4.5898E+01
 EPSSHRINKVR(%)  7.0730E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1603.8302645249994     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -868.67943796126121     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1603.830       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.59E-01  1.00E+00  5.85E-01  1.03E+00  7.68E-01  1.03E+00  1.46E+00  1.47E-01  9.56E-01  9.12E-01  8.94E-01
 


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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       34.249
Stop Time:
Sun Oct 24 04:23:47 CDT 2021
