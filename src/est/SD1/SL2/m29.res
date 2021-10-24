Sat Oct 23 15:04:51 CDT 2021
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
$DATA ../../../../data/SD1/SL2/dat29.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      997
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      897
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1459.37475797044        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7668E+02  5.0913E+00  7.2790E+01  1.3896E+02  1.8392E+02  1.8311E+01 -1.3748E+02 -1.7018E+02 -1.1014E+02 -3.6828E+01
            -4.3629E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2865.68372159782        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0411E+00  1.1910E+00  9.8186E-01  8.8455E-01  1.0581E+00  1.1130E+00  1.5084E+00  9.8100E-01  1.1538E+00  1.1112E+00
             2.0801E+00
 PARAMETER:  1.4024E-01  2.7479E-01  8.1697E-02 -2.2680E-02  1.5650E-01  2.0703E-01  5.1108E-01  8.0816E-02  2.4310E-01  2.0543E-01
             8.3243E-01
 GRADIENT:   2.2312E+02  3.9687E+01 -2.0541E+01 -9.9113E-01  1.8761E+01  3.4630E+01  4.3677E+01 -1.5874E+00  4.9266E+00  9.5891E+00
            -2.3841E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2884.17664937029        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0005E+00  1.3017E+00  1.0670E+00  8.2243E-01  1.1707E+00  1.0476E+00  1.2192E+00  1.2010E-01  1.0645E+00  1.2172E+00
             2.3140E+00
 PARAMETER:  1.0047E-01  3.6364E-01  1.6483E-01 -9.5490E-02  2.5759E-01  1.4651E-01  2.9819E-01 -2.0194E+00  1.6249E-01  2.9659E-01
             9.3897E-01
 GRADIENT:   8.3877E+01  2.4194E+01 -6.2079E-01 -1.5845E+01  2.0371E+01  7.7615E+00  3.4128E+00 -7.8645E-02 -8.4671E+00  1.3464E-01
             3.0911E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2891.89296392203        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      336
 NPARAMETR:  1.0032E+00  1.6422E+00  1.1633E+00  6.7733E-01  1.4126E+00  1.0727E+00  1.0910E+00  7.1632E-02  1.2499E+00  1.4239E+00
             2.3292E+00
 PARAMETER:  1.0324E-01  5.9606E-01  2.5123E-01 -2.8960E-01  4.4542E-01  1.7020E-01  1.8713E-01 -2.5362E+00  3.2305E-01  4.5337E-01
             9.4551E-01
 GRADIENT:   1.2087E+00  6.8935E+00  1.9570E+00  1.9624E+01 -3.0569E+00  5.2970E-01  1.8720E+00 -2.6249E-02  2.3493E+00 -2.3513E+00
             1.1828E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2896.00461805833        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      512
 NPARAMETR:  1.0023E+00  2.0640E+00  9.6023E-01  3.9939E-01  1.7232E+00  1.0698E+00  9.4689E-01  1.0000E-02  1.3461E+00  1.6731E+00
             2.3113E+00
 PARAMETER:  1.0226E-01  8.2464E-01  5.9422E-02 -8.1782E-01  6.4420E-01  1.6746E-01  4.5423E-02 -4.5323E+00  3.9721E-01  6.1468E-01
             9.3783E-01
 GRADIENT:  -3.7788E-01  9.8720E+00  2.1615E+00  3.9799E+00 -6.7722E-01 -1.1735E+00 -2.7221E+00 -5.1920E-05 -3.0253E+00  7.9142E-04
             2.3834E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2897.53009093476        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      689
 NPARAMETR:  1.0026E+00  2.3171E+00  4.9691E-01  2.2520E-01  1.8372E+00  1.0757E+00  8.6722E-01  1.0000E-02  2.0939E+00  1.7419E+00
             2.2998E+00
 PARAMETER:  1.0264E-01  9.4030E-01 -5.9935E-01 -1.3908E+00  7.0826E-01  1.7296E-01 -4.2461E-02 -7.7070E+00  8.3902E-01  6.5497E-01
             9.3280E-01
 GRADIENT:   4.9511E-01  1.7694E+00  2.1452E+00 -4.2422E-01 -2.4558E+00  4.7463E-01 -9.0990E-01  0.0000E+00  3.5840E-01  9.6707E-01
            -3.7720E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2897.92811132213        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      869
 NPARAMETR:  1.0013E+00  2.4630E+00  2.7253E-01  1.4865E-01  1.8991E+00  1.0716E+00  8.2976E-01  1.0000E-02  2.9222E+00  1.7507E+00
             2.2931E+00
 PARAMETER:  1.0125E-01  1.0014E+00 -1.2000E+00 -1.8062E+00  7.4137E-01  1.6915E-01 -8.6616E-02 -1.0329E+01  1.1723E+00  6.6002E-01
             9.2990E-01
 GRADIENT:  -2.7126E+00  3.7695E+01  1.4828E+00  4.3408E+00  1.1970E+00 -1.4548E+00 -1.3396E+00  0.0000E+00  1.7349E+00 -3.5036E+00
             4.4822E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2898.47055296695        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1050             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0025E+00  2.4631E+00  2.1481E-01  1.2204E-01  1.9127E+00  1.0758E+00  8.2412E-01  1.0000E-02  3.0537E+00  1.7883E+00
             2.2858E+00
 PARAMETER:  1.0246E-01  1.0014E+00 -1.4380E+00 -2.0034E+00  7.4853E-01  1.7306E-01 -9.3442E-02 -1.1471E+01  1.2164E+00  6.8128E-01
             9.2670E-01
 GRADIENT:   9.1416E+01  6.0603E+02  8.6373E-01  8.8127E+00  4.4117E+01  1.9090E+01  4.2321E+00  0.0000E+00  3.5961E+00  1.3494E+01
             1.6209E+01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -2898.47863077216        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1143
 NPARAMETR:  1.0024E+00  2.4650E+00  2.1174E-01  1.2317E-01  1.9152E+00  1.0756E+00  8.2453E-01  1.0000E-02  3.0700E+00  1.7882E+00
             2.2859E+00
 PARAMETER:  1.0241E-01  1.0022E+00 -1.4524E+00 -1.9942E+00  7.4983E-01  1.7292E-01 -9.2942E-02 -1.1471E+01  1.2217E+00  6.8120E-01
             9.2676E-01
 GRADIENT:   1.2572E-01 -3.5674E+00 -2.6560E-02  2.7343E-01  2.2940E-01  7.0013E-02  6.8301E-02  0.0000E+00 -3.7904E-01  1.8338E-01
            -2.0635E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1143
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2572E-03 -1.9956E-02 -2.2848E-05  3.1231E-02 -1.9061E-02
 SE:             2.9520E-02  2.7759E-02  2.1511E-05  1.5324E-02  2.6759E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6603E-01  4.7219E-01  2.8817E-01  4.1546E-02  4.7627E-01

 ETASHRINKSD(%)  1.1043E+00  7.0048E+00  9.9928E+01  4.8662E+01  1.0354E+01
 ETASHRINKVR(%)  2.1965E+00  1.3519E+01  1.0000E+02  7.3645E+01  1.9636E+01
 EBVSHRINKSD(%)  1.1668E+00  6.5100E+00  9.9922E+01  5.9822E+01  6.7028E+00
 EBVSHRINKVR(%)  2.3199E+00  1.2596E+01  1.0000E+02  8.3858E+01  1.2956E+01
 RELATIVEINF(%)  9.7660E+01  3.3315E+01  3.5693E-05  4.8617E+00  6.2041E+01
 EPSSHRINKSD(%)  1.7198E+01
 EPSSHRINKVR(%)  3.1437E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          897
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1648.5757285691827     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2898.4786307721583     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1249.9029022029756     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.66
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2898.479       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  2.46E+00  2.12E-01  1.23E-01  1.92E+00  1.08E+00  8.25E-01  1.00E-02  3.07E+00  1.79E+00  2.29E+00
 


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
 #CPUT: Total CPU Time in Seconds,       77.562
Stop Time:
Sat Oct 23 15:05:05 CDT 2021
