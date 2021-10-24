Sun Oct 24 04:34:37 CDT 2021
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
$DATA ../../../../data/SD4/D2/dat37.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m37.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1538.75632740477        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3623E+02 -7.3131E+01 -2.4600E+01 -3.5904E+01  5.8819E+01  5.8435E+00 -1.2686E+02 -4.4447E+00 -1.0929E+02 -3.0406E+01
             1.8268E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1589.70080021274        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      171
 NPARAMETR:  1.0829E+00  1.1482E+00  1.0593E+00  8.9548E-01  1.1449E+00  1.3203E+00  2.5191E+00  1.0337E+00  1.8092E+00  1.1447E+00
             8.9396E-01
 PARAMETER:  1.7966E-01  2.3819E-01  1.5760E-01 -1.0398E-02  2.3534E-01  3.7789E-01  1.0239E+00  1.3315E-01  6.9286E-01  2.3515E-01
            -1.2092E-02
 GRADIENT:   4.3480E+01 -1.7816E+01 -2.6528E+01 -1.3924E+01  2.6672E+01  5.1579E+01  5.5750E+01  3.3186E+00  5.3957E+01  1.1791E+01
            -1.3440E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1601.15781920056        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  1.0226E+00  6.5653E-01  2.4876E+00  1.4242E+00  1.3026E+00  1.2082E+00  3.2572E+00  2.3999E+00  1.2210E+00  1.6681E+00
             8.7470E-01
 PARAMETER:  1.2231E-01 -3.2078E-01  1.0113E+00  4.5361E-01  3.6437E-01  2.8912E-01  1.2809E+00  9.7541E-01  2.9964E-01  6.1166E-01
            -3.3872E-02
 GRADIENT:  -3.1398E+01  2.3704E+01 -5.3182E+00  3.8149E+01 -2.8599E+00  2.5015E+01  2.3452E+01  1.0522E+01  2.6111E+01  2.9610E+01
            -2.3391E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1617.49216896830        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      524
 NPARAMETR:  1.0672E+00  8.2276E-01  1.2336E+00  1.1777E+00  1.0187E+00  1.1170E+00  2.3563E+00  1.2090E+00  1.0416E+00  1.0351E+00
             9.3835E-01
 PARAMETER:  1.6503E-01 -9.5089E-02  3.0996E-01  2.6355E-01  1.1849E-01  2.1066E-01  9.5710E-01  2.8983E-01  1.4080E-01  1.3452E-01
             3.6369E-02
 GRADIENT:   3.3090E+01 -2.9461E+00 -6.9765E+00  7.2263E+00  6.5992E+00 -7.2869E+00 -7.0044E+00  3.9275E+00 -1.2153E+00  2.0824E+00
             4.5361E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1618.34064180701        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      699
 NPARAMETR:  1.0486E+00  8.1780E-01  1.1057E+00  1.1663E+00  9.5531E-01  1.1366E+00  2.4611E+00  9.9955E-01  1.0372E+00  9.5049E-01
             9.2848E-01
 PARAMETER:  1.4747E-01 -1.0114E-01  2.0045E-01  2.5388E-01  5.4283E-02  2.2800E-01  1.0006E+00  9.9552E-02  1.3652E-01  4.9223E-02
             2.5798E-02
 GRADIENT:   4.5451E-01  4.3030E-01  1.7632E-01  2.5255E-01 -2.1227E-01 -1.3336E-02  1.0118E-02 -6.9323E-02  1.3516E-01 -3.0128E-02
            -7.5351E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1618.35097008948        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      878             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0527E+00  8.0455E-01  1.1139E+00  1.1703E+00  9.5515E-01  1.1407E+00  2.5157E+00  1.0099E+00  1.0341E+00  9.5244E-01
             9.2860E-01
 PARAMETER:  1.5133E-01 -1.1747E-01  2.0784E-01  2.5724E-01  5.4112E-02  2.3167E-01  1.0226E+00  1.0981E-01  1.3353E-01  5.1277E-02
             2.5918E-02
 GRADIENT:   7.8878E+02  3.4243E+01  5.4332E+00  3.4018E+02  8.2405E+00  1.6061E+02  1.7857E+02  5.5279E-01  2.0458E+01  1.1010E+00
             9.0122E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1618.35743433111        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1035
 NPARAMETR:  1.0509E+00  8.0622E-01  1.1137E+00  1.1715E+00  9.5512E-01  1.1392E+00  2.5036E+00  1.0053E+00  1.0328E+00  9.5051E-01
             9.2847E-01
 PARAMETER:  1.4962E-01 -1.1540E-01  2.0772E-01  2.5826E-01  5.4079E-02  2.3035E-01  1.0177E+00  1.0524E-01  1.3223E-01  4.9243E-02
             2.5786E-02
 GRADIENT:   4.5477E+00  3.4612E-01  4.7755E-01 -2.2336E+00  6.4531E-02  9.6252E-01  1.7483E+00 -1.4131E-01  6.8182E-03 -1.3079E-01
            -7.1695E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1618.35796988183        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1205
 NPARAMETR:  1.0509E+00  8.0534E-01  1.1133E+00  1.1714E+00  9.5478E-01  1.1392E+00  2.5041E+00  1.0058E+00  1.0327E+00  9.5107E-01
             9.2852E-01
 PARAMETER:  1.4962E-01 -1.1649E-01  2.0729E-01  2.5822E-01  5.3725E-02  2.3036E-01  1.0179E+00  1.0583E-01  1.3218E-01  4.9837E-02
             2.5836E-02
 GRADIENT:   3.2392E-03  1.2246E-01  3.1864E-01 -1.1802E-01 -6.2876E-04  3.0881E-05 -2.1216E-02 -2.1674E-02 -1.2160E-02  1.2103E-02
            -1.4333E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1205
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.5416E-04  2.0243E-02 -4.2968E-02 -2.4324E-02 -2.0710E-02
 SE:             2.9922E-02  2.3073E-02  1.4446E-02  2.2104E-02  2.0051E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7989E-01  3.8030E-01  2.9363E-03  2.7114E-01  3.0166E-01

 ETASHRINKSD(%)  1.0000E-10  2.2704E+01  5.1604E+01  2.5948E+01  3.2827E+01
 ETASHRINKVR(%)  1.0000E-10  4.0253E+01  7.6578E+01  4.5164E+01  5.4877E+01
 EBVSHRINKSD(%)  3.0408E-01  2.2372E+01  5.5899E+01  2.5282E+01  2.9625E+01
 EBVSHRINKVR(%)  6.0723E-01  3.9739E+01  8.0551E+01  4.4172E+01  5.0473E+01
 RELATIVEINF(%)  9.9006E+01  1.1681E+01  4.2851E+00  1.0456E+01  1.1503E+01
 EPSSHRINKSD(%)  4.5892E+01
 EPSSHRINKVR(%)  7.0723E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1618.3579698818298     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -883.20714331809165     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1618.358       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  8.05E-01  1.11E+00  1.17E+00  9.55E-01  1.14E+00  2.50E+00  1.01E+00  1.03E+00  9.51E-01  9.29E-01
 


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
 #CPUT: Total CPU Time in Seconds,       41.823
Stop Time:
Sun Oct 24 04:34:46 CDT 2021
