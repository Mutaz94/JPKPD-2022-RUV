Sat Sep 18 11:42:59 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat48.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
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
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m48.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1640.13591125151        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6808E+01 -3.9764E+01  2.2025E+00 -4.9513E+01  4.1341E+00 -2.8672E+01 -5.7780E+00 -3.3122E+00  7.0398E+00  8.7305E+00
            -4.2774E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1645.59981105282        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8872E-01  1.0837E+00  9.0921E-01  9.7010E-01  9.8777E-01  1.0630E+00  1.1465E+00  1.0958E+00  9.1822E-01  8.4640E-01
             1.1109E+00
 PARAMETER:  8.8661E-02  1.8041E-01  4.8196E-03  6.9645E-02  8.7692E-02  1.6110E-01  2.3672E-01  1.9152E-01  1.4677E-02 -6.6767E-02
             2.0521E-01
 GRADIENT:   3.1275E+00 -9.8264E+00 -6.9530E+00 -1.0472E+01  1.4510E+01 -1.2641E+00  5.1476E+00  2.1460E+00 -1.1686E+00  2.4343E-01
             4.5715E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1646.15858531590        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0012E+00  1.0986E+00  8.1169E-01  9.6553E-01  9.2712E-01  1.0789E+00  1.1288E+00  9.8005E-01  9.3459E-01  7.7198E-01
             1.0968E+00
 PARAMETER:  1.0118E-01  1.9399E-01 -1.0864E-01  6.4925E-02  2.4333E-02  1.7597E-01  2.2117E-01  7.9844E-02  3.2358E-02 -1.5879E-01
             1.9239E-01
 GRADIENT:   2.7284E+01  5.8670E+00 -2.2209E+00  5.5068E+00  4.2284E-01  5.1698E+00  3.8373E+00  2.0660E+00  1.3692E+00 -1.4027E+00
            -9.6235E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1646.43841635275        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      335
 NPARAMETR:  1.0076E+00  1.1917E+00  7.6752E-01  9.0832E-01  9.5273E-01  1.0837E+00  1.0357E+00  8.7096E-01  9.7909E-01  8.2136E-01
             1.1012E+00
 PARAMETER:  1.0755E-01  2.7540E-01 -1.6459E-01  3.8381E-03  5.1574E-02  1.8040E-01  1.3506E-01 -3.8157E-02  7.8868E-02 -9.6790E-02
             1.9637E-01
 GRADIENT:   1.6061E+00  1.3067E+00  6.7826E-01  1.8724E+00 -1.3689E+00  2.4193E-01 -5.1889E-01 -1.0140E-01 -1.2140E-01  2.9918E-01
             2.2599E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1646.56124006787        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      511
 NPARAMETR:  1.0081E+00  1.4229E+00  5.5508E-01  7.4863E-01  9.5390E-01  1.0844E+00  9.2182E-01  5.6517E-01  1.0954E+00  7.9306E-01
             1.1011E+00
 PARAMETER:  1.0807E-01  4.5267E-01 -4.8864E-01 -1.8951E-01  5.2807E-02  1.8103E-01  1.8595E-02 -4.7063E-01  1.9111E-01 -1.3186E-01
             1.9629E-01
 GRADIENT:   5.0892E-02  2.6549E+00  1.5096E-01  1.3203E+00 -2.7458E+00 -9.5142E-02  1.6345E-02  1.9760E-01  2.9634E-01  4.5153E-01
             1.8415E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1646.57679554647        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      688
 NPARAMETR:  1.0079E+00  1.4828E+00  5.1765E-01  7.0652E-01  9.6948E-01  1.0845E+00  8.9664E-01  4.7766E-01  1.1355E+00  7.9752E-01
             1.1013E+00
 PARAMETER:  1.0789E-01  4.9392E-01 -5.5845E-01 -2.4741E-01  6.9008E-02  1.8114E-01 -9.1003E-03 -6.3886E-01  2.2705E-01 -1.2625E-01
             1.9646E-01
 GRADIENT:  -3.8060E-01 -9.2692E-01 -4.3616E-01 -4.6327E-01  3.3052E-01 -8.5516E-02  2.1770E-01  1.4615E-01  1.8806E-01  1.3365E-01
             6.1602E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1646.61799202869        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      865
 NPARAMETR:  1.0081E+00  1.4942E+00  4.7975E-01  6.9536E-01  9.5016E-01  1.0852E+00  8.9430E-01  1.7533E-01  1.1386E+00  7.8306E-01
             1.1009E+00
 PARAMETER:  1.0810E-01  5.0160E-01 -6.3450E-01 -2.6333E-01  4.8872E-02  1.8174E-01 -1.1709E-02 -1.6411E+00  2.2982E-01 -1.4455E-01
             1.9616E-01
 GRADIENT:  -9.9026E-02 -2.0297E-01 -8.1385E-01  6.0982E-01 -7.5964E-02  7.3519E-02  3.3148E-01  2.8541E-02  3.0553E-01  4.1303E-01
             1.9430E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1646.62888408728        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1040
 NPARAMETR:  1.0082E+00  1.4993E+00  4.7754E-01  6.9184E-01  9.5224E-01  1.0850E+00  8.9238E-01  3.9587E-02  1.1427E+00  7.8450E-01
             1.1010E+00
 PARAMETER:  1.0815E-01  5.0498E-01 -6.3910E-01 -2.6840E-01  5.1060E-02  1.8155E-01 -1.3864E-02 -3.1293E+00  2.3340E-01 -1.4271E-01
             1.9624E-01
 GRADIENT:  -7.5369E-03  6.7092E-03 -1.8168E-01  1.2931E-01 -1.0120E-01  1.5703E-03  1.5532E-01  1.0337E-03  1.1739E-01  1.2859E-01
             6.4751E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1646.62953498946        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1215
 NPARAMETR:  1.0082E+00  1.4992E+00  4.7848E-01  6.9196E-01  9.5296E-01  1.0850E+00  8.9200E-01  1.0000E-02  1.1426E+00  7.8492E-01
             1.1010E+00
 PARAMETER:  1.0815E-01  5.0491E-01 -6.3715E-01 -2.6823E-01  5.1822E-02  1.8154E-01 -1.4288E-02 -4.6541E+00  2.3331E-01 -1.4217E-01
             1.9625E-01
 GRADIENT:   1.1599E-03  9.1630E-04  3.7276E-03 -4.7810E-03 -3.7010E-03 -6.7088E-04  7.1114E-04  0.0000E+00 -1.6717E-04 -8.7163E-05
            -2.8600E-05

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1646.62953498946        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1237
 NPARAMETR:  1.0082E+00  1.4992E+00  4.7848E-01  6.9196E-01  9.5296E-01  1.0850E+00  8.9200E-01  1.0000E-02  1.1426E+00  7.8492E-01
             1.1010E+00
 PARAMETER:  1.0815E-01  5.0491E-01 -6.3715E-01 -2.6823E-01  5.1822E-02  1.8154E-01 -1.4288E-02 -4.6541E+00  2.3331E-01 -1.4217E-01
             1.9625E-01
 GRADIENT:   1.1599E-03  9.1630E-04  3.7276E-03 -4.7810E-03 -3.7010E-03 -6.7088E-04  7.1114E-04  0.0000E+00 -1.6717E-04 -8.7163E-05
            -2.8600E-05

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1237
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9182E-04 -1.8157E-02 -3.1930E-04  1.4891E-02 -2.7010E-02
 SE:             2.9832E-02  2.4463E-02  1.3467E-04  2.3643E-02  2.0462E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9487E-01  4.5794E-01  1.7742E-02  5.2880E-01  1.8684E-01

 ETASHRINKSD(%)  5.9647E-02  1.8046E+01  9.9549E+01  2.0792E+01  3.1449E+01
 ETASHRINKVR(%)  1.1926E-01  3.2836E+01  9.9998E+01  3.7260E+01  5.3008E+01
 EBVSHRINKSD(%)  4.4626E-01  1.7908E+01  9.9604E+01  2.1397E+01  3.1285E+01
 EBVSHRINKVR(%)  8.9053E-01  3.2609E+01  9.9998E+01  3.8216E+01  5.2783E+01
 RELATIVEINF(%)  9.9088E+01  4.3544E+00  1.3001E-04  3.8642E+00  5.7388E+00
 EPSSHRINKSD(%)  4.3767E+01
 EPSSHRINKVR(%)  6.8378E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1646.6295349894563     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -911.47870842571808     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.66
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.94
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1646.630       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.50E+00  4.78E-01  6.92E-01  9.53E-01  1.08E+00  8.92E-01  1.00E-02  1.14E+00  7.85E-01  1.10E+00
 


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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.23E+02
 
 TH 2
+       -5.77E+00  4.32E+02
 
 TH 3
+        4.70E+00  2.61E+02  8.24E+02
 
 TH 4
+       -1.46E+01  2.88E+02 -5.41E+02  1.11E+03
 
 TH 5
+       -4.74E+00 -3.24E+02 -7.29E+02  4.55E+02  9.42E+02
 
 TH 6
+        1.25E+00 -1.00E+00  2.38E+00 -4.92E+00 -7.20E-01  1.67E+02
 
 TH 7
+       -5.44E-01  1.72E+01 -3.95E+01 -8.99E+00  2.52E+00  6.92E-01  1.21E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.48E+00 -1.95E+01 -4.30E+01  5.60E+01 -8.46E+00 -2.34E-01  1.59E+01  0.00E+00  6.67E+01
 
 TH10
+       -1.10E+00 -1.46E+01 -4.71E+01 -1.71E+01 -7.44E+01 -1.10E+00  2.03E+01  0.00E+00  1.21E+01  8.71E+01
 
 TH11
+       -6.89E+00 -1.51E+01 -2.94E+01 -1.65E+00 -6.78E+00  2.27E+00  9.33E+00  0.00E+00  7.98E+00  1.98E+01  1.76E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       20.655
Stop Time:
Sat Sep 18 11:43:21 CDT 2021
