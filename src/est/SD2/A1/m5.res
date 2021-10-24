Sat Oct 23 17:19:20 CDT 2021
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
$DATA ../../../../data/SD2/A1/dat5.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m5.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2393.67504650606        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1690E+02  5.1403E+01  3.9234E+01  4.4060E+01  3.1261E+01  3.4187E+01 -2.2190E+01  1.1960E+01  2.1319E+01 -1.6315E+01
            -1.0227E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2590.64352224795        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.8970E-01  9.4990E-01  9.0677E-01  1.0851E+00  9.2930E-01  9.5661E-01  1.1802E+00  7.1904E-01  8.2017E-01  9.7018E-01
             1.4616E+00
 PARAMETER:  8.9650E-02  4.8603E-02  2.1383E-03  1.8166E-01  2.6671E-02  5.5640E-02  2.6570E-01 -2.2984E-01 -9.8243E-02  6.9729E-02
             4.7955E-01
 GRADIENT:   2.8288E+02  5.3359E+01 -1.6197E+01  1.4556E+02  3.4114E+01  7.3089E+00 -1.3863E+00  1.0052E+01 -1.3689E+01 -3.5560E-01
            -7.2693E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2595.62203451407        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.8573E-01  8.0654E-01  8.2796E-01  1.1275E+00  8.2991E-01  9.6287E-01  1.1322E+00  3.8223E-01  8.8395E-01  9.6276E-01
             1.4610E+00
 PARAMETER:  8.5632E-02 -1.1500E-01 -8.8785E-02  2.1996E-01 -8.6444E-02  6.2160E-02  2.2420E-01 -8.6173E-01 -2.3356E-02  6.2044E-02
             4.7912E-01
 GRADIENT:   2.7197E+02 -5.3959E+00 -1.9545E+01  1.2098E+02  5.7049E+01  1.0882E+01 -9.5104E+00  3.8162E+00  3.9070E+00  8.5287E+00
            -6.4858E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2599.18915001269        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  9.2888E-01  8.5234E-01  8.3322E-01  1.0809E+00  8.4842E-01  9.3139E-01  1.1514E+00  2.7080E-01  8.9568E-01  9.6293E-01
             1.4857E+00
 PARAMETER:  2.6220E-02 -5.9772E-02 -8.2456E-02  1.7779E-01 -6.4376E-02  2.8923E-02  2.4094E-01 -1.2064E+00 -1.0177E-02  6.2222E-02
             4.9586E-01
 GRADIENT:   1.2191E+02 -2.0653E+00 -5.6046E+00  4.7835E+01  2.9917E+01  1.9476E+00 -4.2879E+00  1.7837E+00  4.8343E+00  6.7715E+00
            -3.5250E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2602.44963900826        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      382
 NPARAMETR:  9.4996E-01  7.6329E-01  7.6361E-01  1.1398E+00  7.5079E-01  9.6545E-01  1.3658E+00  1.6361E-01  8.6103E-01  8.3272E-01
             1.5311E+00
 PARAMETER:  4.8663E-02 -1.7012E-01 -1.6970E-01  2.3085E-01 -1.8663E-01  6.4839E-02  4.1176E-01 -1.7102E+00 -4.9628E-02 -8.3052E-02
             5.2598E-01
 GRADIENT:   2.6273E+00  5.0767E+00  1.7285E+00  4.7445E+00 -3.1895E+00  1.8639E+00  6.7362E-01  5.9808E-01 -4.7099E-01 -6.9869E-01
             8.2476E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2602.80246207851        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      557
 NPARAMETR:  9.4879E-01  7.4594E-01  7.5021E-01  1.1433E+00  7.3672E-01  9.6075E-01  1.3682E+00  3.5869E-02  8.6352E-01  8.2571E-01
             1.5220E+00
 PARAMETER:  4.7435E-02 -1.9311E-01 -1.8740E-01  2.3390E-01 -2.0554E-01  5.9956E-02  4.1350E-01 -3.2279E+00 -4.6743E-02 -9.1507E-02
             5.2000E-01
 GRADIENT:   5.5120E-03  5.0797E-01  2.5434E+00 -1.2227E+00 -2.4568E+00 -3.6088E-02 -5.4778E-01  2.9435E-02  1.4627E-01 -4.9850E-01
            -1.8922E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2602.82328121446        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      732
 NPARAMETR:  9.4879E-01  7.3870E-01  7.4128E-01  1.1462E+00  7.2922E-01  9.6090E-01  1.3803E+00  1.0000E-02  8.6246E-01  8.1908E-01
             1.5217E+00
 PARAMETER:  4.7436E-02 -2.0286E-01 -1.9937E-01  2.3645E-01 -2.1577E-01  6.0112E-02  4.2229E-01 -4.7242E+00 -4.7971E-02 -9.9577E-02
             5.1983E-01
 GRADIENT:   6.9569E-03  1.4405E-02  2.9695E-02 -3.4325E-02 -2.5267E-02  2.6987E-03 -1.2027E-02  0.0000E+00 -1.9486E-02 -4.5673E-02
             3.4361E-02

0ITERATION NO.:   31    OBJECTIVE VALUE:  -2602.82328121446        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      754
 NPARAMETR:  9.4879E-01  7.3870E-01  7.4128E-01  1.1462E+00  7.2922E-01  9.6090E-01  1.3803E+00  1.0000E-02  8.6246E-01  8.1908E-01
             1.5217E+00
 PARAMETER:  4.7436E-02 -2.0286E-01 -1.9937E-01  2.3645E-01 -2.1577E-01  6.0112E-02  4.2229E-01 -4.7242E+00 -4.7971E-02 -9.9577E-02
             5.1983E-01
 GRADIENT:   6.9569E-03  1.4405E-02  2.9695E-02 -3.4325E-02 -2.5267E-02  2.6987E-03 -1.2027E-02  0.0000E+00 -1.9486E-02 -4.5673E-02
             3.4361E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      754
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.3278E-04  2.1344E-03 -3.4386E-04 -3.6991E-03 -2.3674E-03
 SE:             2.9723E-02  2.1768E-02  1.8320E-04  2.6765E-02  2.3324E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7496E-01  9.2189E-01  6.0518E-02  8.9008E-01  9.1915E-01

 ETASHRINKSD(%)  4.2279E-01  2.7073E+01  9.9386E+01  1.0334E+01  2.1862E+01
 ETASHRINKVR(%)  8.4379E-01  4.6817E+01  9.9996E+01  1.9601E+01  3.8945E+01
 EBVSHRINKSD(%)  7.0520E-01  2.6411E+01  9.9335E+01  1.0774E+01  2.2493E+01
 EBVSHRINKVR(%)  1.4054E+00  4.5846E+01  9.9996E+01  2.0386E+01  3.9927E+01
 RELATIVEINF(%)  9.8586E+01  1.4009E+01  9.8133E-04  3.5933E+01  6.8589E+00
 EPSSHRINKSD(%)  2.3191E+01
 EPSSHRINKVR(%)  4.1004E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          700
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1286.5139464865417     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2602.8232812144633     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1316.3093347279216     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2602.823       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.49E-01  7.39E-01  7.41E-01  1.15E+00  7.29E-01  9.61E-01  1.38E+00  1.00E-02  8.62E-01  8.19E-01  1.52E+00
 


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
 #CPUT: Total CPU Time in Seconds,       42.903
Stop Time:
Sat Oct 23 17:19:30 CDT 2021
