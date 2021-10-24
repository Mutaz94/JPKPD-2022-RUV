Sat Oct 23 17:44:56 CDT 2021
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
$DATA ../../../../data/SD2/A2/dat56.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m56.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2011.48647538572        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9934E+02  1.2552E+01  1.2263E+02 -1.6792E+01  8.4004E+01  3.7095E+01 -4.3604E+01 -6.7122E+01 -1.2956E+00 -2.0908E+01
            -1.7375E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2481.45007177517        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.6306E-01  9.8423E-01  6.8793E-01  1.0582E+00  8.2332E-01  1.0716E+00  1.2507E+00  8.4750E-01  9.0798E-01  6.1684E-01
             1.8777E+00
 PARAMETER:  6.2359E-02  8.4108E-02 -2.7407E-01  1.5659E-01 -9.4412E-02  1.6912E-01  3.2369E-01 -6.5461E-02  3.4706E-03 -3.8314E-01
             7.3005E-01
 GRADIENT:   2.6302E+01  1.6562E+01 -1.1972E+01  3.4313E+01  2.8678E+01  4.3001E+01  7.9096E+00  1.1213E+01  5.4454E+00  1.2324E+00
             1.3948E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2489.42812059577        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      204
 NPARAMETR:  1.0183E+00  6.8582E-01  4.9523E-01  1.1560E+00  5.4029E-01  9.6571E-01  1.4828E+00  5.8782E-01  8.2293E-01  3.9264E-01
             1.8571E+00
 PARAMETER:  1.1813E-01 -2.7714E-01 -6.0273E-01  2.4496E-01 -5.1565E-01  6.5107E-02  4.9394E-01 -4.3133E-01 -9.4880E-02 -8.3486E-01
             7.1903E-01
 GRADIENT:   2.6383E+01  4.5085E+01  4.6968E+01 -8.3032E+01 -8.7232E+01 -1.4401E+01 -1.5525E+00  4.2704E+00 -1.0569E+01 -8.1228E+00
             1.1836E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2496.54788705655        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      382
 NPARAMETR:  1.0098E+00  6.4582E-01  4.7698E-01  1.1968E+00  5.3515E-01  1.0358E+00  1.4063E+00  3.5066E-01  8.1761E-01  6.6220E-01
             1.8272E+00
 PARAMETER:  1.0972E-01 -3.3724E-01 -6.4028E-01  2.7968E-01 -5.2520E-01  1.3522E-01  4.4093E-01 -9.4795E-01 -1.0137E-01 -3.1218E-01
             7.0277E-01
 GRADIENT:   6.4214E+00  8.0066E+00 -3.4172E+01 -4.4897E+00  3.6464E+01  1.2565E+01  2.2768E+00  3.2825E+00 -1.4441E+01  8.8012E+00
             8.5945E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2506.90679362732        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      561
 NPARAMETR:  9.9740E-01  3.5855E-01  2.7081E-01  1.2044E+00  3.0199E-01  9.9640E-01  1.1596E+00  5.1913E-02  9.3494E-01  7.7743E-01
             1.7557E+00
 PARAMETER:  9.7396E-02 -9.2568E-01 -1.2063E+00  2.8599E-01 -1.0973E+00  9.6393E-02  2.4806E-01 -2.8582E+00  3.2732E-02 -1.5177E-01
             6.6287E-01
 GRADIENT:  -2.2531E+01 -6.6815E+00  6.5594E+00  3.6593E+01 -6.5930E+00 -3.8099E+00  4.9155E+00  3.7292E-02  5.6067E+00  2.7860E+00
             9.8441E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2507.89690958280        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      738
 NPARAMETR:  1.0072E+00  3.4786E-01  2.5258E-01  1.1693E+00  2.8811E-01  1.0055E+00  1.0815E+00  4.8816E-02  9.2708E-01  7.8462E-01
             1.7420E+00
 PARAMETER:  1.0720E-01 -9.5597E-01 -1.2760E+00  2.5643E-01 -1.1444E+00  1.0545E-01  1.7835E-01 -2.9197E+00  2.4280E-02 -1.4255E-01
             6.5505E-01
 GRADIENT:  -7.5676E-01 -9.2800E-01  6.1627E-01 -3.9568E-01 -6.9442E-02  4.6650E-02  5.6412E-01  2.0082E-02  6.4419E-01  1.1154E+00
             6.2747E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2507.90963952458        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      917
 NPARAMETR:  1.0075E+00  3.4862E-01  2.5277E-01  1.1695E+00  2.8833E-01  1.0053E+00  1.0730E+00  1.9913E-02  9.2539E-01  7.8255E-01
             1.7423E+00
 PARAMETER:  1.0749E-01 -9.5377E-01 -1.2753E+00  2.5662E-01 -1.1436E+00  1.0532E-01  1.7043E-01 -3.8164E+00  2.2462E-02 -1.4519E-01
             6.5519E-01
 GRADIENT:  -1.4272E-01 -4.3717E-02  1.4387E-01 -2.2474E-01 -6.0415E-01 -3.1131E-03  3.6421E-02  3.3580E-03  5.6942E-02  8.7299E-02
             6.6463E-02

0ITERATION NO.:   34    OBJECTIVE VALUE:  -2507.91072541028        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1055
 NPARAMETR:  1.0077E+00  3.4860E-01  2.5259E-01  1.1696E+00  2.8835E-01  1.0054E+00  1.0707E+00  1.0000E-02  9.2515E-01  7.8231E-01
             1.7423E+00
 PARAMETER:  1.0770E-01 -9.5267E-01 -1.2746E+00  2.5662E-01 -1.1439E+00  1.0537E-01  1.6997E-01 -5.3728E+00  2.2219E-02 -1.4517E-01
             6.5525E-01
 GRADIENT:  -2.8648E-03  5.8067E-01  1.2035E+00 -2.6325E-02 -6.2382E-01  9.9749E-04  3.4818E-02  0.0000E+00  3.4602E-03  4.9949E-02
             2.4394E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1055
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.0595E-04  6.4361E-03 -7.8265E-05 -4.5530E-03  2.5391E-03
 SE:             2.9680E-02  1.8823E-02  2.4606E-04  2.8526E-02  2.6417E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8102E-01  7.3240E-01  7.5043E-01  8.7319E-01  9.2343E-01

 ETASHRINKSD(%)  5.6832E-01  3.6941E+01  9.9176E+01  4.4358E+00  1.1500E+01
 ETASHRINKVR(%)  1.1334E+00  6.0236E+01  9.9993E+01  8.6748E+00  2.1678E+01
 EBVSHRINKSD(%)  8.1231E-01  3.6236E+01  9.9168E+01  4.1420E+00  1.1909E+01
 EBVSHRINKVR(%)  1.6180E+00  5.9342E+01  9.9993E+01  8.1124E+00  2.2400E+01
 RELATIVEINF(%)  9.8369E+01  9.1445E+00  6.3876E-04  5.7775E+01  5.7128E+00
 EPSSHRINKSD(%)  2.4476E+01
 EPSSHRINKVR(%)  4.2962E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          700
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1286.5139464865417     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2507.9107254102832     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1221.3967789237415     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.86
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2507.911       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  3.49E-01  2.53E-01  1.17E+00  2.88E-01  1.01E+00  1.07E+00  1.00E-02  9.25E-01  7.83E-01  1.74E+00
 


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
 #CPUT: Total CPU Time in Seconds,       63.090
Stop Time:
Sat Oct 23 17:45:09 CDT 2021
