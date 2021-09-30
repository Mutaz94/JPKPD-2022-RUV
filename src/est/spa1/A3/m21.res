Thu Sep 30 00:04:04 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat21.csv ignore=@
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
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   9.29837857795722        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5391E+02  3.7047E+01  2.3594E+02 -5.3195E+01  1.7342E+02  5.5106E+01 -6.3211E+01 -1.9165E+02 -3.7088E+01 -1.1987E+02
            -3.8273E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1490.63335284949        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0837E+00  1.0341E+00  8.8188E-01  1.2062E+00  9.5847E-01  8.4426E-01  9.8038E-01  1.0057E+00  9.5071E-01  1.0009E+00
             5.2460E+00
 PARAMETER:  1.8035E-01  1.3354E-01 -2.5695E-02  2.8746E-01  5.7578E-02 -6.9293E-02  8.0184E-02  1.0568E-01  4.9456E-02  1.0086E-01
             1.7575E+00
 GRADIENT:   5.7861E+00  1.5313E+01 -1.9926E+01  5.4436E+01 -8.8425E+00 -1.0374E+01  1.1284E+01  8.6432E+00  2.3557E+01  2.2499E+01
             3.4025E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1539.17108053520        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0640E+00  6.9619E-01  3.2586E-01  1.3146E+00  4.0258E-01  9.2378E-01  5.1639E-01  9.5917E-02  1.2379E+00  2.1814E-01
             4.4824E+00
 PARAMETER:  1.6202E-01 -2.6214E-01 -1.0213E+00  3.7352E-01 -8.0986E-01  2.0715E-02 -5.6090E-01 -2.2443E+00  3.1342E-01 -1.4226E+00
             1.6002E+00
 GRADIENT:  -2.0934E+01  9.1898E+01  5.7742E+01  1.3283E+02 -1.1009E+02  8.4522E-01 -4.8199E-01  2.3087E-02  2.8174E+01 -1.4026E-01
             2.5645E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1599.95241515164        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      276
 NPARAMETR:  1.0248E+00  5.7224E-01  1.8340E-01  1.1650E+00  2.9002E-01  9.4904E-01  1.6403E-01  1.8728E-02  1.4302E+00  7.1613E-01
             3.0726E+00
 PARAMETER:  1.2453E-01 -4.5819E-01 -1.5961E+00  2.5269E-01 -1.1378E+00  4.7692E-02 -1.7077E+00 -3.8777E+00  4.5784E-01 -2.3389E-01
             1.2225E+00
 GRADIENT:  -6.3434E+01  5.4643E+01 -7.8037E+01  7.9375E+01  4.9891E+01 -4.0694E+00  3.5779E-01 -1.0345E-03  9.3227E+00  1.5013E+01
             6.2506E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1618.49621890865        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      452
 NPARAMETR:  1.0420E+00  4.6813E-01  1.9434E-01  1.0996E+00  2.6753E-01  9.4866E-01  1.2449E-01  3.3756E-02  1.3328E+00  7.5015E-01
             2.6841E+00
 PARAMETER:  1.4113E-01 -6.5902E-01 -1.5382E+00  1.9490E-01 -1.2185E+00  4.7291E-02 -1.9836E+00 -3.2886E+00  3.8730E-01 -1.8749E-01
             1.0874E+00
 GRADIENT:  -1.4097E+01  3.4196E+01  7.9103E+00  3.8960E+00 -6.1504E+00 -1.3038E+00  1.8042E-01 -8.2813E-03  3.3807E+00  1.8317E+01
            -3.6572E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1625.15196687734        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      628
 NPARAMETR:  1.0488E+00  3.5820E-01  1.5961E-01  1.0602E+00  2.2249E-01  9.4852E-01  1.5430E-02  7.2608E-02  1.3618E+00  6.7043E-01
             2.7657E+00
 PARAMETER:  1.4764E-01 -9.2666E-01 -1.7350E+00  1.5845E-01 -1.4029E+00  4.7151E-02 -4.0714E+00 -2.5227E+00  4.0882E-01 -2.9984E-01
             1.1173E+00
 GRADIENT:  -5.2826E-01 -1.0025E-01 -6.9308E-01  8.0667E-01  1.8222E+00 -2.8752E-01  2.6222E-03 -1.1605E-01  1.7711E-01 -4.7818E-01
             7.8093E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1630.72610851354        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      810
 NPARAMETR:  1.0467E+00  3.4460E-01  1.4198E-01  1.0322E+00  2.0921E-01  9.5034E-01  1.0000E-02  9.0194E-01  1.4389E+00  5.3775E-01
             2.6320E+00
 PARAMETER:  1.4562E-01 -9.6537E-01 -1.8521E+00  1.3166E-01 -1.4644E+00  4.9061E-02 -4.9332E+00 -3.2103E-03  4.6385E-01 -5.2036E-01
             1.0677E+00
 GRADIENT:   4.7684E-01  2.7891E+00 -1.1544E+01 -2.4171E+00  6.0996E+00 -5.8848E-01  0.0000E+00 -6.4912E+00 -6.5721E-01 -6.2021E-01
            -6.5896E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1632.00581109899        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      988
 NPARAMETR:  1.0467E+00  3.5293E-01  1.5610E-01  1.0627E+00  2.2053E-01  9.5178E-01  1.0000E-02  1.1530E+00  1.4255E+00  4.1178E-01
             2.6278E+00
 PARAMETER:  1.4560E-01 -9.4149E-01 -1.7573E+00  1.6080E-01 -1.4117E+00  5.0583E-02 -4.6878E+00  2.4238E-01  4.5453E-01 -7.8728E-01
             1.0662E+00
 GRADIENT:  -5.7147E-02  4.5977E-01  7.4148E-01  1.2823E-01 -1.3777E+00  5.1647E-02  0.0000E+00 -6.9780E-01  4.6263E-02 -2.2615E-01
            -7.3931E-01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1632.00921018434        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1045
 NPARAMETR:  1.0467E+00  3.5254E-01  1.5576E-01  1.0621E+00  2.2033E-01  9.5169E-01  1.0000E-02  1.1627E+00  1.4266E+00  4.1188E-01
             2.6274E+00
 PARAMETER:  1.4560E-01 -9.4259E-01 -1.7594E+00  1.6024E-01 -1.4126E+00  5.0489E-02 -4.7030E+00  2.5075E-01  4.5528E-01 -7.8701E-01
             1.0660E+00
 GRADIENT:   1.3844E-02 -4.2182E-02  1.2896E-01 -6.2235E-02  6.2818E-05  3.5265E-02  0.0000E+00  7.0969E-03  2.9840E-02  3.1446E-02
            -1.8495E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1045
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7419E-03 -2.4377E-04  1.7051E-02 -4.9591E-03  2.0398E-02
 SE:             2.9148E-02  1.4663E-04  2.1031E-02  2.7116E-02  1.5607E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5235E-01  9.6406E-02  4.1750E-01  8.5489E-01  1.9120E-01

 ETASHRINKSD(%)  2.3518E+00  9.9509E+01  2.9544E+01  9.1590E+00  4.7716E+01
 ETASHRINKVR(%)  4.6482E+00  9.9998E+01  5.0360E+01  1.7479E+01  7.2664E+01
 EBVSHRINKSD(%)  2.3680E+00  9.9501E+01  2.8913E+01  7.2459E+00  4.8750E+01
 EBVSHRINKVR(%)  4.6799E+00  9.9998E+01  4.9466E+01  1.3967E+01  7.3734E+01
 RELATIVEINF(%)  9.5167E+01  4.3135E-04  6.3859E+00  5.1453E+01  1.4670E+00
 EPSSHRINKSD(%)  2.8489E+01
 EPSSHRINKVR(%)  4.8862E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1632.0092101843366     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -713.07067697966386     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.03
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1632.009       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  3.53E-01  1.56E-01  1.06E+00  2.20E-01  9.52E-01  1.00E-02  1.16E+00  1.43E+00  4.12E-01  2.63E+00
 


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
+        1.07E+03
 
 TH 2
+       -9.31E+01  8.13E+04
 
 TH 3
+       -1.97E+05 -8.61E+04  1.68E+04
 
 TH 4
+       -9.85E+00  9.17E+02  7.08E+02  4.16E+02
 
 TH 5
+        1.78E+02 -7.43E+03 -2.05E+04 -6.83E+02  3.74E+04
 
 TH 6
+        2.99E+00 -1.08E+02 -9.46E+01 -5.99E+00  2.79E+01  1.98E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        4.07E-01 -8.82E+04  2.40E+03 -7.44E-01  4.26E+01  1.08E+00  0.00E+00  3.35E+01
 
 TH 9
+        8.64E+00 -3.99E+04  5.54E+02 -5.43E+00  1.96E+02  6.57E-01  0.00E+00 -2.02E+00  6.73E+01
 
 TH10
+       -1.66E-02 -6.07E+02 -9.27E+04  1.44E-01  5.55E+02  3.20E+00  0.00E+00  3.62E+01  1.52E+01  9.48E+01
 
 TH11
+       -1.78E+01 -9.53E+03  1.04E+04 -3.03E+00 -9.75E+03  2.21E+00  0.00E+00  1.17E+01  4.89E+00  1.03E+01  6.26E+01
 
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
 #CPUT: Total CPU Time in Seconds,       24.121
Stop Time:
Thu Sep 30 00:04:32 CDT 2021
