Sat Oct 23 17:35:43 CDT 2021
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
$DATA ../../../../data/SD2/A2/dat3.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m3.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2205.00166289715        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9406E+02  2.2147E+01  9.7956E+01 -5.3001E+00  1.1356E+02  4.8062E+01 -3.0435E+01 -1.2021E+02 -5.0002E+00 -2.2546E+00
            -1.3839E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2558.97695665003        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0291E+00  1.0173E+00  8.1688E-01  1.0582E+00  8.6627E-01  1.0630E+00  9.1713E-01  1.0961E+00  1.0190E+00  6.2684E-01
             1.7112E+00
 PARAMETER:  1.2871E-01  1.1718E-01 -1.0226E-01  1.5660E-01 -4.3554E-02  1.6111E-01  1.3490E-02  1.9173E-01  1.1884E-01 -3.6706E-01
             6.3722E-01
 GRADIENT:   2.0240E+02  5.3990E+01  1.1785E+01  5.8494E+01 -6.6341E+00  5.2441E+01 -6.2104E+00 -9.8437E+00  2.5674E+00  1.0078E+00
            -3.8276E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2560.70155479805        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0302E+00  8.3276E-01  6.6861E-01  1.1671E+00  6.8739E-01  1.0134E+00  1.2127E+00  1.0772E+00  9.3773E-01  5.3091E-01
             1.6597E+00
 PARAMETER:  1.2980E-01 -8.3015E-02 -3.0255E-01  2.5455E-01 -2.7486E-01  1.1333E-01  2.9288E-01  1.7440E-01  3.5704E-02 -5.3316E-01
             6.0661E-01
 GRADIENT:   2.1466E+02  7.3917E+01  1.1422E+01  1.4216E+02  7.1486E+00  2.8098E+01  1.3288E+01  7.2047E+00 -9.7854E-01  6.8297E+00
            -5.8765E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2564.71569326411        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      335
 NPARAMETR:  1.0266E+00  7.9187E-01  7.0325E-01  1.1779E+00  6.9290E-01  9.8993E-01  1.2376E+00  1.0193E+00  9.1480E-01  4.8685E-01
             1.7332E+00
 PARAMETER:  1.2628E-01 -1.3336E-01 -2.5205E-01  2.6375E-01 -2.6687E-01  8.9874E-02  3.1317E-01  1.1912E-01  1.0950E-02 -6.1980E-01
             6.4998E-01
 GRADIENT:   1.1649E+01  3.0225E+01  7.5994E+00  3.0008E+01 -1.1899E+01  3.5796E+00  7.3859E+00  3.1187E+00 -6.8974E+00  3.3079E+00
            -8.3397E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2575.97222653463        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      511
 NPARAMETR:  9.9301E-01  4.1877E-01  3.4960E-01  1.2086E+00  3.7016E-01  9.8820E-01  1.4354E+00  7.8659E-01  9.4358E-01  2.7146E-01
             1.6843E+00
 PARAMETER:  9.2987E-02 -7.7043E-01 -9.5098E-01  2.8945E-01 -8.9381E-01  8.8127E-02  4.6146E-01 -1.4004E-01  4.1924E-02 -1.2039E+00
             6.2136E-01
 GRADIENT:  -6.2290E+01 -2.5520E+01  1.6554E+00  2.3645E+00  1.9129E+01 -2.2958E+00 -9.6352E+00 -1.1166E+01 -1.2828E+01 -3.1670E+00
             5.3672E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2581.86998006781        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      687
 NPARAMETR:  1.0197E+00  3.8583E-01  2.9004E-01  1.1765E+00  3.2454E-01  9.9306E-01  1.5081E+00  9.6868E-01  1.0537E+00  2.4944E-01
             1.5456E+00
 PARAMETER:  1.1947E-01 -8.5237E-01 -1.1377E+00  2.6257E-01 -1.0253E+00  9.3040E-02  5.1087E-01  6.8182E-02  1.5226E-01 -1.2885E+00
             5.3540E-01
 GRADIENT:   7.4404E-01  2.4366E+00 -2.8930E+00  1.3103E+00  3.0326E+00  5.4320E-01  2.0615E-01 -1.2990E+00  9.0061E-01 -1.0855E+00
            -3.0419E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2582.09799122555        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      863
 NPARAMETR:  1.0193E+00  3.7643E-01  2.8394E-01  1.1729E+00  3.1795E-01  9.9170E-01  1.4975E+00  9.4964E-01  1.0526E+00  3.1517E-01
             1.5448E+00
 PARAMETER:  1.1910E-01 -8.7701E-01 -1.1590E+00  2.5950E-01 -1.0459E+00  9.1665E-02  5.0377E-01  4.8328E-02  1.5126E-01 -1.0546E+00
             5.3487E-01
 GRADIENT:  -1.2721E-01 -2.6037E-01 -1.0794E-01 -3.0467E-01  4.7931E-01  2.7218E-02  1.1704E-01  2.1907E-02  6.5457E-02 -2.6574E-02
             6.7461E-02

0ITERATION NO.:   31    OBJECTIVE VALUE:  -2582.09799122555        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  1.0193E+00  3.7643E-01  2.8394E-01  1.1729E+00  3.1795E-01  9.9170E-01  1.4975E+00  9.4964E-01  1.0526E+00  3.1517E-01
             1.5448E+00
 PARAMETER:  1.1910E-01 -8.7701E-01 -1.1590E+00  2.5950E-01 -1.0459E+00  9.1665E-02  5.0377E-01  4.8328E-02  1.5126E-01 -1.0546E+00
             5.3487E-01
 GRADIENT:  -1.2721E-01 -2.6037E-01 -1.0794E-01 -3.0467E-01  4.7931E-01  2.7218E-02  1.1704E-01  2.1907E-02  6.5457E-02 -2.6574E-02
             6.7461E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      885
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.9472E-04  1.6609E-02 -1.6644E-02 -5.7771E-03  3.2787E-03
 SE:             2.9761E-02  2.4398E-02  2.2547E-02  2.9040E-02  1.2102E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7602E-01  4.9603E-01  4.6040E-01  8.4231E-01  7.8646E-01

 ETASHRINKSD(%)  2.9758E-01  1.8264E+01  2.4465E+01  2.7140E+00  5.9456E+01
 ETASHRINKVR(%)  5.9427E-01  3.3193E+01  4.2945E+01  5.3544E+00  8.3562E+01
 EBVSHRINKSD(%)  6.7231E-01  1.6508E+01  2.4735E+01  3.3629E+00  6.0344E+01
 EBVSHRINKVR(%)  1.3401E+00  3.0291E+01  4.3352E+01  6.6127E+00  8.4274E+01
 RELATIVEINF(%)  9.8648E+01  2.5601E+01  9.6828E+00  8.1597E+01  1.8068E+00
 EPSSHRINKSD(%)  2.7071E+01
 EPSSHRINKVR(%)  4.6814E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          700
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1286.5139464865417     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2582.0979912255521     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1295.5840447390103     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.28
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2582.098       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  3.76E-01  2.84E-01  1.17E+00  3.18E-01  9.92E-01  1.50E+00  9.50E-01  1.05E+00  3.15E-01  1.54E+00
 


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
 #CPUT: Total CPU Time in Seconds,       53.273
Stop Time:
Sat Oct 23 17:35:54 CDT 2021
