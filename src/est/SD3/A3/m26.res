Sat Oct 23 22:32:46 CDT 2021
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
$DATA ../../../../data/SD3/A3/dat26.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   419.805939301306        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4751E+02  9.8973E+01  2.5224E+02  3.2144E+00  3.0189E+02  6.9925E+01 -1.0536E+02 -4.0749E+02 -2.0501E+02 -1.8586E+02
            -4.0889E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1448.14641584815        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0160E+00  1.0036E+00  9.4842E-01  1.0733E+00  9.2392E-01  7.0334E-01  1.0122E+00  1.0746E+00  1.0196E+00  9.2911E-01
             4.7691E+00
 PARAMETER:  1.1587E-01  1.0355E-01  4.7044E-02  1.7070E-01  2.0867E-02 -2.5191E-01  1.1208E-01  1.7191E-01  1.1942E-01  2.6476E-02
             1.6622E+00
 GRADIENT:  -1.7544E+01 -8.0620E+00 -1.5380E+01  1.0944E+00 -8.3780E+00 -2.8067E+01  7.0976E+00  8.9927E+00  1.4864E+01  2.3173E+01
             2.3323E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1465.31012703933        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0050E+00  5.5814E-01  4.6452E-01  1.3786E+00  4.4344E-01  7.6860E-01  1.2939E+00  6.5163E-01  9.8149E-01  3.3705E-01
             4.3353E+00
 PARAMETER:  1.0500E-01 -4.8315E-01 -6.6674E-01  4.2108E-01 -7.1319E-01 -1.6319E-01  3.5769E-01 -3.2828E-01  8.1321E-02 -9.8751E-01
             1.5668E+00
 GRADIENT:  -6.6986E+01  6.6860E+01  2.0055E+01  2.1159E+02 -5.8892E+01 -1.8155E+01 -7.2045E-01  9.4810E+00 -1.2826E+01  4.8272E+00
             1.7556E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1491.28727797830        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.8108E-01  5.1719E-01  2.0672E-01  1.2296E+00  2.8103E-01  9.4101E-01  1.0539E+00  9.3789E-01  1.0252E+00  1.1264E-01
             3.0778E+00
 PARAMETER:  8.0897E-02 -5.5935E-01 -1.4764E+00  3.0672E-01 -1.1693E+00  3.9194E-02  1.5250E-01  3.5876E-02  1.2485E-01 -2.0836E+00
             1.2242E+00
 GRADIENT:  -4.6503E+01  9.5029E+01  1.0133E+00  2.6889E+02 -2.5576E+01  2.6966E+01 -1.4622E+01 -4.8896E+00 -1.1062E+02 -5.9957E-01
            -2.2359E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1511.30683531637        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      411
 NPARAMETR:  9.9764E-01  4.5652E-01  2.2056E-01  1.1887E+00  2.7525E-01  8.6080E-01  1.3392E+00  1.1247E+00  1.1700E+00  4.2020E-02
             2.7123E+00
 PARAMETER:  9.7634E-02 -6.8413E-01 -1.4116E+00  2.7290E-01 -1.1901E+00 -4.9891E-02  3.9211E-01  2.1748E-01  2.5699E-01 -3.0696E+00
             1.0978E+00
 GRADIENT:  -1.7215E+01  4.3822E+01  2.4223E+01  1.3294E+02 -4.6207E+01  1.0289E+00  2.4723E-01 -2.8665E+00 -6.0533E+01 -5.8589E-02
            -1.0995E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1532.92761706012        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      586
 NPARAMETR:  1.0026E+00  3.8810E-01  1.7054E-01  1.0485E+00  2.3350E-01  8.4881E-01  1.0993E+00  1.1457E+00  1.5517E+00  1.0000E-02
             2.9641E+00
 PARAMETER:  1.0258E-01 -8.4650E-01 -1.6688E+00  1.4737E-01 -1.3546E+00 -6.3922E-02  1.9470E-01  2.3598E-01  5.3932E-01 -5.0731E+00
             1.1866E+00
 GRADIENT:  -1.4123E+00  7.5835E+00  3.3528E+00  6.5971E+00 -1.1217E+01 -5.7708E-02 -4.2883E-01  2.2752E+00  2.4263E+00  0.0000E+00
            -1.3432E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1533.09211789168        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      762
 NPARAMETR:  1.0033E+00  3.7715E-01  1.7025E-01  1.0400E+00  2.3210E-01  8.4933E-01  1.1457E+00  1.1110E+00  1.5246E+00  1.0000E-02
             2.9707E+00
 PARAMETER:  1.0326E-01 -8.7512E-01 -1.6705E+00  1.3925E-01 -1.3606E+00 -6.3313E-02  2.3602E-01  2.0524E-01  5.2174E-01 -5.1722E+00
             1.1888E+00
 GRADIENT:   2.8197E-01 -1.1335E-01 -4.0315E-01 -1.3542E-01  9.4355E-01  1.5152E-03 -3.8354E-02 -1.2478E-01  8.5401E-02  0.0000E+00
            -2.9201E-01

0ITERATION NO.:   31    OBJECTIVE VALUE:  -1533.09211789168        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      784
 NPARAMETR:  1.0033E+00  3.7715E-01  1.7025E-01  1.0400E+00  2.3210E-01  8.4933E-01  1.1457E+00  1.1110E+00  1.5246E+00  1.0000E-02
             2.9707E+00
 PARAMETER:  1.0326E-01 -8.7512E-01 -1.6705E+00  1.3925E-01 -1.3606E+00 -6.3313E-02  2.3602E-01  2.0524E-01  5.2174E-01 -5.1722E+00
             1.1888E+00
 GRADIENT:   2.8197E-01 -1.1335E-01 -4.0315E-01 -1.3542E-01  9.4355E-01  1.5152E-03 -3.8354E-02 -1.2478E-01  8.5401E-02  0.0000E+00
            -2.9201E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      784
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3485E-03  1.5813E-03  6.9610E-04 -8.5225E-03  5.3129E-04
 SE:             2.8846E-02  1.5832E-02  2.0340E-02  2.7134E-02  4.5491E-04
 N:                     100         100         100         100         100

 P VAL.:         9.6271E-01  9.2044E-01  9.7270E-01  7.5346E-01  2.4284E-01

 ETASHRINKSD(%)  3.3621E+00  4.6959E+01  3.1859E+01  9.0970E+00  9.8476E+01
 ETASHRINKVR(%)  6.6112E+00  7.1867E+01  5.3568E+01  1.7366E+01  9.9977E+01
 EBVSHRINKSD(%)  3.5823E+00  4.7406E+01  3.1341E+01  7.4022E+00  9.8773E+01
 EBVSHRINKVR(%)  7.0363E+00  7.2339E+01  5.2859E+01  1.4257E+01  9.9985E+01
 RELATIVEINF(%)  9.2747E+01  6.0088E+00  7.6980E+00  6.1791E+01  1.1482E-03
 EPSSHRINKSD(%)  2.7567E+01
 EPSSHRINKVR(%)  4.7534E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1533.0921178916772     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -614.15358468700447     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.79
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1533.092       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  3.77E-01  1.70E-01  1.04E+00  2.32E-01  8.49E-01  1.15E+00  1.11E+00  1.52E+00  1.00E-02  2.97E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       66.048
Stop Time:
Sat Oct 23 22:32:58 CDT 2021
