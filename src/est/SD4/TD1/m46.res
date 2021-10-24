Sun Oct 24 03:47:25 CDT 2021
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
$DATA ../../../../data/SD4/TD1/dat46.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m46.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1666.62412630941        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5802E+02 -7.7352E+01 -4.4142E+01 -4.4431E+01  3.9770E+01  4.8228E+01  1.7153E+00  9.4406E+00  2.5626E+01  3.0375E+01
            -2.2274E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1675.71246933689        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.7820E-01  1.1723E+00  1.2504E+00  9.7441E-01  1.1171E+00  9.3716E-01  9.7625E-01  9.4894E-01  8.0283E-01  7.5420E-01
             1.1710E+00
 PARAMETER:  7.7959E-02  2.5896E-01  3.2346E-01  7.4073E-02  2.1078E-01  3.5103E-02  7.5960E-02  4.7588E-02 -1.1961E-01 -1.8210E-01
             2.5782E-01
 GRADIENT:   2.3592E+01  3.2068E+00  2.7731E+01 -2.5829E+01  3.3683E+00 -1.0542E+01 -8.2686E+00 -1.4759E+01 -1.6214E+01 -1.8733E+01
             1.6891E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1676.81058764512        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.6380E-01  9.4553E-01  1.3092E+00  1.1237E+00  1.0224E+00  9.5735E-01  1.1090E+00  1.2594E+00  7.3291E-01  6.2010E-01
             1.1729E+00
 PARAMETER:  6.3126E-02  4.3994E-02  3.6943E-01  2.1664E-01  1.2217E-01  5.6413E-02  2.0345E-01  3.3064E-01 -2.1073E-01 -3.7788E-01
             2.5948E-01
 GRADIENT:  -1.0484E+01  3.3973E+00  1.4429E+01 -9.4763E+00 -8.6999E+00 -1.5738E+00 -1.1251E+01 -4.3121E-01 -1.5767E+01 -1.6714E+01
             1.9283E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1682.13466910779        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  9.6798E-01  8.7003E-01  1.4784E+00  1.1874E+00  1.0789E+00  9.6303E-01  1.1997E+00  1.2458E+00  7.9702E-01  9.0171E-01
             1.0635E+00
 PARAMETER:  6.7458E-02 -3.9231E-02  4.9097E-01  2.7178E-01  1.7593E-01  6.2332E-02  2.8204E-01  3.1975E-01 -1.2687E-01 -3.4605E-03
             1.6153E-01
 GRADIENT:   4.8190E+00  7.1687E+00  1.5731E+00  7.6799E+00 -5.1566E+00  5.5628E-01  8.6134E-01 -1.6495E-01  5.5885E-01  1.8823E+00
            -3.8536E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1682.33957923036        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  9.6150E-01  7.0505E-01  1.7309E+00  1.3016E+00  1.0912E+00  9.6065E-01  1.3304E+00  1.3793E+00  7.5478E-01  9.1251E-01
             1.0830E+00
 PARAMETER:  6.0740E-02 -2.4949E-01  6.4862E-01  3.6357E-01  1.8727E-01  5.9850E-02  3.8544E-01  4.2158E-01 -1.8133E-01  8.4410E-03
             1.7975E-01
 GRADIENT:  -7.6023E+00  8.4753E+00  3.9864E+00  1.2712E+01 -7.5625E+00  3.2552E-01  5.0314E-02 -8.1628E-01 -8.7078E-01  1.4753E-01
             2.7991E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1682.64386622733        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  9.6276E-01  5.0519E-01  1.9189E+00  1.4344E+00  1.0817E+00  9.5725E-01  1.5535E+00  1.4777E+00  7.2126E-01  9.2765E-01
             1.0703E+00
 PARAMETER:  6.2050E-02 -5.8282E-01  7.5176E-01  4.6072E-01  1.7856E-01  5.6310E-02  5.4053E-01  4.9048E-01 -2.2675E-01  2.4903E-02
             1.6797E-01
 GRADIENT:   4.4293E-01  5.6943E+00  7.4484E-01  1.6022E+01 -3.6332E+00 -3.5853E-01 -4.6459E-01 -2.5970E-01 -1.2124E+00  5.6019E-01
            -1.5902E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1682.75058328448        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1055
 NPARAMETR:  9.6190E-01  3.7400E-01  2.0959E+00  1.5207E+00  1.0864E+00  9.5658E-01  1.8031E+00  1.5919E+00  7.0178E-01  9.3598E-01
             1.0703E+00
 PARAMETER:  6.1151E-02 -8.8350E-01  8.3999E-01  5.1918E-01  1.8289E-01  5.5610E-02  6.8949E-01  5.6492E-01 -2.5414E-01  3.3834E-02
             1.6797E-01
 GRADIENT:   1.4692E+00  3.2496E+00 -3.6789E-01  1.2004E+01 -1.1032E+00 -1.3344E-01 -3.3955E-01 -3.9003E-03 -3.8360E-01  2.7656E-01
            -1.4164E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1682.80977652075        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1238             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6139E-01  3.3326E-01  2.1606E+00  1.5410E+00  1.0879E+00  9.5656E-01  1.9666E+00  1.6314E+00  6.9432E-01  9.3695E-01
             1.0731E+00
 PARAMETER:  6.0623E-02 -9.9883E-01  8.7039E-01  5.3243E-01  1.8425E-01  5.5587E-02  7.7629E-01  5.8942E-01 -2.6482E-01  3.4870E-02
             1.7058E-01
 GRADIENT:   3.3734E+02  3.3830E+01  7.1132E+00  6.3499E+02  6.7401E+00  3.4999E+01  6.5627E+00  1.9824E+00  1.8251E+01  7.7585E-01
             1.2229E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1682.81244031593        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1420             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6134E-01  3.3130E-01  2.1577E+00  1.5411E+00  1.0881E+00  9.5655E-01  1.9801E+00  1.6323E+00  6.9125E-01  9.3288E-01
             1.0734E+00
 PARAMETER:  6.0572E-02 -1.0047E+00  8.6906E-01  5.3248E-01  1.8441E-01  5.5575E-02  7.8317E-01  5.8998E-01 -2.6926E-01  3.0523E-02
             1.7087E-01
 GRADIENT:   3.3703E+02  3.3186E+01  6.4543E+00  6.3304E+02  8.7315E+00  3.4987E+01  6.5588E+00  2.0692E+00  1.7549E+01  2.6690E-01
             1.1988E+00

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1682.81330618979        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:     1549
 NPARAMETR:  9.6135E-01  3.3132E-01  2.1579E+00  1.5409E+00  1.0878E+00  9.5656E-01  1.9911E+00  1.6321E+00  6.9203E-01  9.3397E-01
             1.0736E+00
 PARAMETER:  6.0582E-02 -1.0047E+00  8.6916E-01  5.3236E-01  1.8412E-01  5.5585E-02  7.8870E-01  5.8990E-01 -2.6813E-01  3.1688E-02
             1.7106E-01
 GRADIENT:   1.2312E+00  1.9776E-01  9.3657E-02 -8.5099E+00  1.0180E-01  9.6811E-02  6.7203E-02  6.5639E-02  1.3057E-01  4.1822E-02
            -1.1404E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1549
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3754E-04  3.5371E-03 -3.3967E-02 -8.7615E-03 -4.3597E-02
 SE:             2.9806E-02  1.1910E-02  1.7375E-02  2.6755E-02  1.9477E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9632E-01  7.6648E-01  5.0581E-02  7.4331E-01  2.5195E-02

 ETASHRINKSD(%)  1.4649E-01  6.0099E+01  4.1793E+01  1.0368E+01  3.4750E+01
 ETASHRINKVR(%)  2.9276E-01  8.4079E+01  6.6120E+01  1.9662E+01  5.7424E+01
 EBVSHRINKSD(%)  5.2372E-01  6.2097E+01  4.6197E+01  1.0177E+01  3.1700E+01
 EBVSHRINKVR(%)  1.0447E+00  8.5634E+01  7.1053E+01  1.9318E+01  5.3351E+01
 RELATIVEINF(%)  9.6915E+01  3.3393E-01  7.0495E+00  2.0639E+00  8.7866E+00
 EPSSHRINKSD(%)  4.3317E+01
 EPSSHRINKVR(%)  6.7871E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1682.8133061897950     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -947.66247962605678     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.44
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1682.813       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.61E-01  3.31E-01  2.16E+00  1.54E+00  1.09E+00  9.57E-01  1.99E+00  1.63E+00  6.92E-01  9.34E-01  1.07E+00
 


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
 #CPUT: Total CPU Time in Seconds,       47.489
Stop Time:
Sun Oct 24 03:47:35 CDT 2021
