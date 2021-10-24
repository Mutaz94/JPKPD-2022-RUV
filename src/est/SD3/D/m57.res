Sun Oct 24 01:01:15 CDT 2021
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
$DATA ../../../../data/SD3/D/dat57.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m57.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2131.39687739431        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5141E+02 -4.9920E+01 -5.3676E+01 -1.3766E+01  6.1869E+01  4.3316E+01 -6.2425E+00  1.6984E+01 -7.6847E+00  3.1197E+01
            -2.8850E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2145.57408571790        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      192
 NPARAMETR:  1.0382E+00  1.0375E+00  1.2004E+00  1.0518E+00  1.0451E+00  1.0237E+00  1.0441E+00  8.9618E-01  1.0704E+00  7.9691E-01
             1.0692E+00
 PARAMETER:  1.3746E-01  1.3680E-01  2.8265E-01  1.5055E-01  1.4413E-01  1.2344E-01  1.4311E-01 -9.6152E-03  1.6802E-01 -1.2702E-01
             1.6690E-01
 GRADIENT:   2.6086E-01 -1.5755E+00  1.0753E+00 -3.9474E+00  2.0064E+01  4.2130E+00 -2.9248E+00 -1.0449E-01  3.9056E+00 -4.5048E+00
             1.3098E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2148.00859823787        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0527E+00  8.1191E-01  1.0543E+00  1.2038E+00  8.7742E-01  9.8696E-01  1.4395E+00  7.1434E-01  9.5737E-01  6.6434E-01
             1.0557E+00
 PARAMETER:  1.5132E-01 -1.0837E-01  1.5290E-01  2.8546E-01 -3.0767E-02  8.6878E-02  4.6430E-01 -2.3640E-01  5.6438E-02 -3.0896E-01
             1.5420E-01
 GRADIENT:   3.5360E+01  1.8385E+01  5.5608E+00  3.0820E+01 -7.0765E+00 -1.0725E+01  1.3216E+00  6.5190E-01  1.2048E+01 -5.6885E+00
             7.0751E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2149.86080911527        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.0291E+00  6.9856E-01  1.0229E+00  1.2504E+00  8.2493E-01  1.0181E+00  1.5958E+00  6.1715E-01  8.5814E-01  7.0000E-01
             1.0386E+00
 PARAMETER:  1.2872E-01 -2.5873E-01  1.2268E-01  3.2347E-01 -9.2461E-02  1.1794E-01  5.6737E-01 -3.8264E-01 -5.2988E-02 -2.5667E-01
             1.3789E-01
 GRADIENT:  -1.2739E+01  1.0538E+01  1.0417E+01  4.8601E+00 -1.9575E+01  2.6750E+00 -8.4285E-01 -3.5059E-01 -3.1200E+00 -6.4091E-01
            -1.5092E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2151.01679149509        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0322E+00  3.8226E-01  1.2653E+00  1.4650E+00  8.3581E-01  1.0005E+00  2.1588E+00  8.1812E-01  8.1778E-01  7.7735E-01
             1.0357E+00
 PARAMETER:  1.3166E-01 -8.6166E-01  3.3529E-01  4.8187E-01 -7.9358E-02  1.0051E-01  8.6954E-01 -1.0074E-01 -1.0117E-01 -1.5187E-01
             1.3505E-01
 GRADIENT:   5.9673E+00  6.5289E+00  3.3641E+00  2.0264E+01 -6.5115E+00 -1.3835E+00  3.0542E-01 -8.1551E-01  2.7121E-03  1.8673E-01
            -9.1342E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2151.31839754114        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      902
 NPARAMETR:  1.0274E+00  2.3280E-01  1.3777E+00  1.5610E+00  8.3722E-01  9.9907E-01  2.7557E+00  9.3501E-01  7.9043E-01  7.9112E-01
             1.0359E+00
 PARAMETER:  1.2700E-01 -1.3576E+00  4.2045E-01  5.4533E-01 -7.7665E-02  9.9072E-02  1.1137E+00  3.2801E-02 -1.3518E-01 -1.3431E-01
             1.3526E-01
 GRADIENT:   1.5554E+00  3.9877E+00  2.8775E+00  1.8374E+01 -6.4486E+00 -6.2893E-01 -3.1628E-03 -1.8064E-01 -7.8442E-01 -2.5454E-02
            -1.7492E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2151.44520158407        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1077
 NPARAMETR:  1.0246E+00  1.5065E-01  1.4447E+00  1.6124E+00  8.3899E-01  9.9926E-01  3.3714E+00  1.0033E+00  7.7811E-01  8.0193E-01
             1.0352E+00
 PARAMETER:  1.2428E-01 -1.7928E+00  4.6791E-01  5.7771E-01 -7.5556E-02  9.9263E-02  1.3153E+00  1.0327E-01 -1.5088E-01 -1.2074E-01
             1.3463E-01
 GRADIENT:  -1.2707E+00  2.5596E+00  2.7572E+00  1.2853E+01 -6.8635E+00  1.7512E-01  3.0983E-01  3.8781E-01 -5.8919E-01  2.4612E-01
            -1.6085E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2151.61641769296        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1253
 NPARAMETR:  1.0224E+00  6.7926E-02  1.5008E+00  1.6679E+00  8.3586E-01  9.9786E-01  4.7334E+00  1.0519E+00  7.6204E-01  8.0633E-01
             1.0344E+00
 PARAMETER:  1.2216E-01 -2.5893E+00  5.0603E-01  6.1154E-01 -7.9292E-02  9.7856E-02  1.6546E+00  1.5060E-01 -1.7176E-01 -1.1526E-01
             1.3385E-01
 GRADIENT:  -2.8637E+00 -4.5875E-02  5.2919E+00  2.5049E+01 -1.2430E+01  3.4751E-01 -4.0953E+00  7.5796E-01  2.7030E-01  1.5382E+00
            -4.0979E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2152.34974770340        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1435             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0239E+00  3.5468E-02  1.4814E+00  1.6543E+00  8.2951E-01  9.9743E-01  6.3168E+00  1.0288E+00  7.5485E-01  7.9244E-01
             1.0345E+00
 PARAMETER:  1.2366E-01 -3.2391E+00  4.9299E-01  6.0337E-01 -8.6926E-02  9.7427E-02  1.9432E+00  1.2838E-01 -1.8124E-01 -1.3264E-01
             1.3389E-01
 GRADIENT:   4.9407E+02  8.4976E+00  1.1634E+01  1.1474E+03  1.2681E+01  4.1483E+01  3.6791E+01 -4.3811E-02  2.3319E+01 -3.9131E-01
             1.3529E+00

0ITERATION NO.:   44    OBJECTIVE VALUE:  -2152.41059249303        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     1569
 NPARAMETR:  1.0235E+00  3.5410E-02  1.4799E+00  1.6646E+00  8.2813E-01  9.9661E-01  6.3118E+00  1.0301E+00  7.5441E-01  7.9511E-01
             1.0345E+00
 PARAMETER:  1.2321E-01 -3.2385E+00  4.9203E-01  6.0973E-01 -8.8091E-02  9.6650E-02  1.9437E+00  1.2906E-01 -1.8193E-01 -1.3058E-01
             1.3388E-01
 GRADIENT:  -8.2124E-01  1.4770E+01  2.0284E-01  4.1681E+00  2.7711E+00  1.2301E-01  2.1850E+01 -1.3783E-01 -2.7152E-01 -4.7830E-01
            -1.0303E-01
 NUMSIGDIG:         2.5         2.3         3.2         2.7         1.5         2.5         2.4         1.6         2.4         1.2
                    3.0

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1569
 NO. OF SIG. DIGITS IN FINAL EST.:  1.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.0548E-04  7.4799E-03 -2.1724E-02 -6.5470E-03 -2.7834E-02
 SE:             2.9871E-02  7.0033E-03  1.7219E-02  2.9196E-02  2.0769E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8383E-01  2.8549E-01  2.0709E-01  8.2257E-01  1.8018E-01

 ETASHRINKSD(%)  1.0000E-10  7.6538E+01  4.2313E+01  2.1914E+00  3.0422E+01
 ETASHRINKVR(%)  1.0000E-10  9.4495E+01  6.6722E+01  4.3347E+00  5.1590E+01
 EBVSHRINKSD(%)  3.4273E-01  8.4303E+01  4.3968E+01  2.5693E+00  2.9233E+01
 EBVSHRINKVR(%)  6.8428E-01  9.7536E+01  6.8605E+01  5.0725E+00  4.9920E+01
 RELATIVEINF(%)  9.9030E+01  1.0001E+00  5.8893E+00  4.1608E+01  9.2285E+00
 EPSSHRINKSD(%)  3.3233E+01
 EPSSHRINKVR(%)  5.5422E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2152.4105924930304     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1233.4720592883577     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.90
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2152.411       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  3.55E-02  1.48E+00  1.66E+00  8.29E-01  9.97E-01  6.32E+00  1.03E+00  7.54E-01  7.94E-01  1.03E+00
 


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
 #CPUT: Total CPU Time in Seconds,       52.209
Stop Time:
Sun Oct 24 01:01:26 CDT 2021
