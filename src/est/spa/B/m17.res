Sat Sep 18 08:20:12 CDT 2021
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
$DATA ../../../../data/spa/B/dat17.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1652.80965336836        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.7622E+01 -6.3561E+01  2.2375E+01 -1.1054E+02 -3.3236E+01  4.8453E-01 -1.6371E+01  5.4414E-01 -3.8987E-01 -3.0272E+01
            -3.0308E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1665.50641892703        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0549E+00  1.1903E+00  1.0505E+00  9.7172E-01  1.2167E+00  9.9257E-01  1.2323E+00  9.3278E-01  8.8995E-01  1.4548E+00
             1.0713E+00
 PARAMETER:  1.5343E-01  2.7422E-01  1.4929E-01  7.1316E-02  2.9617E-01  9.2542E-02  3.0886E-01  3.0417E-02 -1.6593E-02  4.7490E-01
             1.6885E-01
 GRADIENT:   9.3599E+01  1.7195E+01 -6.1021E+00  1.0677E+01  9.3420E+00 -9.0649E-01  1.0062E+01 -6.0173E-01 -4.6977E+00  1.1164E+01
             8.2166E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1665.70332983171        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0535E+00  1.2109E+00  1.0971E+00  9.6328E-01  1.2421E+00  9.9711E-01  1.1979E+00  9.3866E-01  9.2671E-01  1.4644E+00
             1.0801E+00
 PARAMETER:  1.5211E-01  2.9134E-01  1.9271E-01  6.2592E-02  3.1680E-01  9.7107E-02  2.8060E-01  3.6694E-02  2.3888E-02  4.8143E-01
             1.7701E-01
 GRADIENT:   3.2986E+01  1.3032E+01  4.7815E+00  4.1692E+00  3.1361E+00 -2.5598E+00  8.2183E+00 -3.1790E+00 -3.4124E+00  8.4080E+00
             2.3020E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1667.15250159288        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      402
 NPARAMETR:  1.0392E+00  1.2564E+00  1.0031E+00  9.2596E-01  1.1951E+00  1.0037E+00  1.0554E+00  9.3271E-01  1.0145E+00  1.3279E+00
             1.0666E+00
 PARAMETER:  1.3843E-01  3.2826E-01  1.0308E-01  2.3077E-02  2.7822E-01  1.0371E-01  1.5391E-01  3.0342E-02  1.1437E-01  3.8357E-01
             1.6449E-01
 GRADIENT:   5.3842E-01  2.2339E+00  1.0234E+00  3.6690E+00  3.2162E-01 -1.5193E-01 -2.4908E-01 -5.5620E-01  2.7221E-01 -7.1013E-01
            -7.5349E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1667.39623559249        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      582
 NPARAMETR:  1.0430E+00  1.4916E+00  8.0607E-01  7.7150E-01  1.2310E+00  1.0055E+00  9.3001E-01  7.4922E-01  1.1431E+00  1.3365E+00
             1.0710E+00
 PARAMETER:  1.4210E-01  4.9983E-01 -1.1559E-01 -1.5942E-01  3.0782E-01  1.0553E-01  2.7438E-02 -1.8872E-01  2.3375E-01  3.9007E-01
             1.6858E-01
 GRADIENT:   4.6361E+00  5.0241E+00  3.0590E-01  4.7049E+00 -2.0340E-01 -2.8572E-01 -8.1014E-01 -2.2242E-01 -1.5551E+00  5.5735E-01
             6.0405E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1667.79607690176        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      758
 NPARAMETR:  1.0417E+00  1.8196E+00  5.6226E-01  5.5301E-01  1.3211E+00  1.0085E+00  8.1295E-01  4.2118E-01  1.4387E+00  1.3563E+00
             1.0711E+00
 PARAMETER:  1.4083E-01  6.9862E-01 -4.7578E-01 -4.9238E-01  3.7847E-01  1.0844E-01 -1.0708E-01 -7.6469E-01  4.6372E-01  4.0476E-01
             1.6867E-01
 GRADIENT:  -7.7185E-01  9.1466E+00  4.0768E-01  4.2239E+00 -8.8172E-02  1.8204E-01 -7.0631E-01  2.6987E-02 -3.8220E-01 -5.2656E-01
            -1.7692E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1667.85741278317        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      933
 NPARAMETR:  1.0421E+00  1.9030E+00  4.9950E-01  4.9091E-01  1.3541E+00  1.0081E+00  7.9145E-01  3.2032E-01  1.5465E+00  1.3725E+00
             1.0710E+00
 PARAMETER:  1.4127E-01  7.4343E-01 -5.9414E-01 -6.1149E-01  4.0316E-01  1.0811E-01 -1.3389E-01 -1.0384E+00  5.3602E-01  4.1667E-01
             1.6857E-01
 GRADIENT:  -1.2363E-02 -6.6729E-01 -2.1517E-01 -1.4952E-01  9.5579E-02 -3.8257E-03  3.5676E-02  6.9103E-02  3.7106E-02  5.8364E-02
             6.0830E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1667.87023646448        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1114
 NPARAMETR:  1.0427E+00  1.9182E+00  4.8292E-01  4.8046E-01  1.3564E+00  1.0082E+00  7.8833E-01  1.9912E-01  1.5614E+00  1.3704E+00
             1.0699E+00
 PARAMETER:  1.4180E-01  7.5137E-01 -6.2791E-01 -6.3301E-01  4.0480E-01  1.0818E-01 -1.3784E-01 -1.5139E+00  5.4557E-01  4.1510E-01
             1.6761E-01
 GRADIENT:   1.0207E+00 -1.0926E+00 -6.5031E-01  1.0691E-01  7.9475E-01 -1.6612E-02  3.8983E-03  2.9919E-02 -1.1855E-01 -6.4226E-02
            -3.2880E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1667.88237161585        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1292
 NPARAMETR:  1.0412E+00  1.9063E+00  4.9220E-01  4.8890E-01  1.3510E+00  1.0079E+00  7.9148E-01  1.0505E-01  1.5486E+00  1.3688E+00
             1.0709E+00
 PARAMETER:  1.4038E-01  7.4516E-01 -6.0887E-01 -6.1561E-01  4.0082E-01  1.0783E-01 -1.3385E-01 -2.1534E+00  5.3733E-01  4.1396E-01
             1.6855E-01
 GRADIENT:  -2.0099E+00  7.0440E-01  8.8120E-02  1.2907E-01  2.4985E-02 -1.4005E-01 -8.2807E-02  5.5043E-03 -1.0399E-01 -1.0538E-01
            -1.9870E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1667.88627823892        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1469
 NPARAMETR:  1.0422E+00  1.9066E+00  4.9056E-01  4.8828E-01  1.3503E+00  1.0082E+00  7.9168E-01  2.1181E-02  1.5498E+00  1.3687E+00
             1.0714E+00
 PARAMETER:  1.4134E-01  7.4532E-01 -6.1220E-01 -6.1686E-01  4.0032E-01  1.0822E-01 -1.3359E-01 -3.7546E+00  5.3814E-01  4.1383E-01
             1.6898E-01
 GRADIENT:   6.7182E-02 -1.0540E-01 -4.1353E-02  8.2098E-03  1.2265E-02  1.1240E-02  1.5294E-02  2.4094E-04  1.5706E-02  1.3822E-02
             2.5590E-02

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1667.88637039481        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1561
 NPARAMETR:  1.0422E+00  1.9067E+00  4.9071E-01  4.8825E-01  1.3505E+00  1.0082E+00  7.9161E-01  1.0000E-02  1.5499E+00  1.3688E+00
             1.0714E+00
 PARAMETER:  1.4130E-01  7.4538E-01 -6.1190E-01 -6.1692E-01  4.0047E-01  1.0818E-01 -1.3369E-01 -4.7243E+00  5.3822E-01  4.1392E-01
             1.6895E-01
 GRADIENT:  -1.6518E-02  2.8465E-02  2.6786E-03  8.4217E-03  6.9960E-03 -2.0437E-03 -3.8191E-03  0.0000E+00 -5.5785E-03 -4.6294E-03
            -7.4352E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1561
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.2902E-05 -3.7244E-02 -2.9747E-04  3.4352E-02 -4.8199E-02
 SE:             2.9831E-02  2.3988E-02  9.4704E-05  2.1924E-02  2.2438E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9832E-01  1.2052E-01  1.6835E-03  1.1715E-01  3.1708E-02

 ETASHRINKSD(%)  6.0852E-02  1.9636E+01  9.9683E+01  2.6551E+01  2.4830E+01
 ETASHRINKVR(%)  1.2167E-01  3.5416E+01  9.9999E+01  4.6053E+01  4.3494E+01
 EBVSHRINKSD(%)  4.8635E-01  1.8447E+01  9.9746E+01  2.9973E+01  2.1408E+01
 EBVSHRINKVR(%)  9.7033E-01  3.3491E+01  9.9999E+01  5.0962E+01  3.8234E+01
 RELATIVEINF(%)  9.8951E+01  5.1544E+00  8.6255E-05  3.5522E+00  1.9685E+01
 EPSSHRINKSD(%)  4.4245E+01
 EPSSHRINKVR(%)  6.8914E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1667.8863703948068     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -932.73554383106864     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.29
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.58
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1667.886       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  1.91E+00  4.91E-01  4.88E-01  1.35E+00  1.01E+00  7.92E-01  1.00E-02  1.55E+00  1.37E+00  1.07E+00
 


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
+        9.99E+02
 
 TH 2
+       -5.40E+00  3.48E+02
 
 TH 3
+        9.19E+00  1.28E+02  3.46E+02
 
 TH 4
+       -1.67E+01  3.09E+02 -3.00E+02  9.83E+02
 
 TH 5
+       -1.81E+00 -8.15E+01 -1.57E+02  1.53E+02  2.05E+02
 
 TH 6
+        1.97E-02 -7.96E-01  2.93E+00 -2.72E+00  4.66E-01  1.95E+02
 
 TH 7
+        4.45E+00  4.64E+00 -1.30E+01 -1.54E+01 -6.37E+00 -2.48E+00  1.53E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.17E+00 -1.67E+01 -3.52E+01  6.32E+01  4.99E-01 -3.99E-01  1.44E+01  0.00E+00  3.12E+01
 
 TH10
+        4.88E-01 -1.04E+01 -2.41E+01  3.83E+00 -3.58E+01  1.19E-01  7.66E+00  0.00E+00  4.05E+00  4.81E+01
 
 TH11
+       -7.33E+00 -2.09E+01 -3.79E+01  7.16E+00 -2.49E+00  1.71E+00  1.28E+01  0.00E+00  5.26E+00  7.72E+00  1.88E+02
 
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
 #CPUT: Total CPU Time in Seconds,       26.925
Stop Time:
Sat Sep 18 08:20:41 CDT 2021
