Sat Sep 18 06:20:44 CDT 2021
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
$DATA ../../../../data/int/TD2/dat87.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3380.78721938081        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3969E+01 -5.1387E+01 -2.4833E+00 -4.7338E+01  1.3109E+02 -6.6192E+00 -2.4999E+01 -2.0854E+02 -3.7548E+01 -2.7827E+01
            -7.8291E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3901.24757939016        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.9172E-01  1.0874E+00  1.1030E+00  9.8923E-01  1.0203E+00  1.0195E+00  9.4894E-01  9.6981E-01  1.0611E+00  9.5413E-01
             1.1679E+00
 PARAMETER:  9.1688E-02  1.8377E-01  1.9801E-01  8.9169E-02  1.2011E-01  1.1928E-01  4.7586E-02  6.9344E-02  1.5930E-01  5.3041E-02
             2.5520E-01
 GRADIENT:  -1.7361E+01  3.6590E+00 -9.5310E+00  4.6320E+00  2.4873E+01 -1.3785E+00  5.0134E+00  6.3508E+00  4.5534E+00 -5.6650E+00
             2.5402E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3903.67357259834        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0014E+00  1.0023E+00  9.7095E-01  1.0277E+00  9.3904E-01  1.0123E+00  9.2316E-01  6.0558E-01  9.8766E-01  9.0904E-01
             1.1453E+00
 PARAMETER:  1.0140E-01  1.0227E-01  7.0524E-02  1.2730E-01  3.7100E-02  1.1227E-01  2.0046E-02 -4.0157E-01  8.7583E-02  4.6284E-03
             2.3564E-01
 GRADIENT:   4.9520E+00 -2.3557E+01 -4.2430E+01  1.1311E+01  5.4928E+01 -3.7854E+00 -3.0365E+00 -3.3621E+00 -2.1237E+01 -6.9020E+00
             2.0465E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3919.34354114537        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  9.9559E-01  1.0388E+00  1.0478E+00  1.0068E+00  9.6463E-01  1.0217E+00  9.2225E-01  8.8282E-01  1.0677E+00  9.3835E-01
             1.0223E+00
 PARAMETER:  9.5583E-02  1.3804E-01  1.4669E-01  1.0676E-01  6.3990E-02  1.2143E-01  1.9065E-02 -2.4635E-02  1.6553E-01  3.6373E-02
             1.2205E-01
 GRADIENT:   3.8873E+00 -7.9469E-02 -8.4069E+00  1.6568E+00  6.3017E+00  1.5347E+00 -3.4388E+00  1.7364E+00  1.9768E+00 -4.5402E+00
             4.2226E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3919.41194416872        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  9.9620E-01  1.0397E+00  1.0497E+00  1.0066E+00  9.6559E-01  1.0220E+00  9.2401E-01  8.8351E-01  1.0670E+00  9.4054E-01
             1.0221E+00
 PARAMETER:  9.6194E-02  1.3892E-01  1.4855E-01  1.0656E-01  6.4980E-02  1.2172E-01  2.0971E-02 -2.3855E-02  1.6485E-01  3.8694E-02
             1.2189E-01
 GRADIENT:  -3.9950E+01 -1.0831E+01 -8.6338E+00 -8.5158E+00  8.8161E-01 -6.5695E+00 -3.7038E+00  1.5391E+00  1.2720E-01 -4.5502E+00
             3.6953E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3919.43513882387        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      538
 NPARAMETR:  9.9634E-01  1.0399E+00  1.0505E+00  1.0065E+00  9.6581E-01  1.0220E+00  9.2504E-01  8.8312E-01  1.0666E+00  9.4154E-01
             1.0220E+00
 PARAMETER:  9.6333E-02  1.3914E-01  1.4925E-01  1.0651E-01  6.5210E-02  1.2174E-01  2.2078E-02 -2.4295E-02  1.6451E-01  3.9760E-02
             1.2180E-01
 GRADIENT:   5.5279E+00  5.1020E-01 -7.2052E+00  1.9551E+00  5.6895E+00  1.7307E+00 -2.9075E+00  1.4849E+00  1.7879E+00 -3.8634E+00
             3.7122E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3919.44743636548        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:      733             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9647E-01  1.0399E+00  1.0507E+00  1.0065E+00  9.6581E-01  1.0223E+00  9.2571E-01  8.8308E-01  1.0666E+00  9.4160E-01
             1.0220E+00
 PARAMETER:  9.6469E-02  1.3914E-01  1.4946E-01  1.0651E-01  6.5210E-02  1.2210E-01  2.2809E-02 -2.4343E-02  1.6449E-01  3.9830E-02
             1.2179E-01
 GRADIENT:   5.8454E+00  5.4354E-01 -7.0541E+00  1.9459E+00  5.5614E+00  1.9086E+00 -2.8118E+00  1.4650E+00  1.7974E+00 -3.8162E+00
             3.7167E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3919.45202206976        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      916
 NPARAMETR:  9.9647E-01  1.0399E+00  1.0507E+00  1.0065E+00  9.6581E-01  1.0224E+00  9.2680E-01  8.8308E-01  1.0666E+00  9.4160E-01
             1.0220E+00
 PARAMETER:  9.6469E-02  1.3914E-01  1.4950E-01  1.0651E-01  6.5210E-02  1.2216E-01  2.3985E-02 -2.4342E-02  1.6448E-01  3.9829E-02
             1.2179E-01
 GRADIENT:  -3.9331E+01 -1.0609E+01 -8.0657E+00 -8.4648E+00  4.5350E-01 -6.3691E+00 -3.2696E+00  1.4158E+00  9.3532E-02 -4.2088E+00
             3.5408E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3919.50839343988        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1086            RESET HESSIAN, TYPE II
 NPARAMETR:  9.9662E-01  1.0398E+00  1.0508E+00  1.0065E+00  9.6581E-01  1.0386E+00  9.2691E-01  8.8296E-01  1.0665E+00  9.4168E-01
             1.0221E+00
 PARAMETER:  9.6611E-02  1.3907E-01  1.4958E-01  1.0651E-01  6.5208E-02  1.3791E-01  2.4099E-02 -2.4474E-02  1.6440E-01  3.9914E-02
             1.2184E-01
 GRADIENT:   7.3877E+00  4.6581E-01 -6.9641E+00  1.8745E+00  5.5146E+00  9.3443E+00 -2.6461E+00  1.4474E+00  1.7938E+00 -3.7310E+00
             3.8949E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3919.51816898619        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1279             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9671E-01  1.0398E+00  1.0512E+00  1.0065E+00  9.6577E-01  1.0383E+00  9.2754E-01  8.8277E-01  1.0664E+00  9.4186E-01
             1.0220E+00
 PARAMETER:  9.6709E-02  1.3907E-01  1.4992E-01  1.0651E-01  6.5169E-02  1.3763E-01  2.4785E-02 -2.4689E-02  1.6432E-01  4.0105E-02
             1.2180E-01
 GRADIENT:   7.5709E+00  5.2600E-01 -6.6735E+00  1.8452E+00  5.2403E+00  9.2170E+00 -2.5550E+00  1.4053E+00  1.7917E+00 -3.6588E+00
             3.8285E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3919.52279954998        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1440
 NPARAMETR:  9.9672E-01  1.0398E+00  1.0512E+00  1.0065E+00  9.6577E-01  1.0383E+00  9.2882E-01  8.8277E-01  1.0664E+00  9.4186E-01
             1.0220E+00
 PARAMETER:  9.6710E-02  1.3907E-01  1.4997E-01  1.0651E-01  6.5169E-02  1.3763E-01  2.6163E-02 -2.4688E-02  1.6432E-01  4.0104E-02
             1.2181E-01
 GRADIENT:   7.5745E+00  5.9841E-01 -6.6664E+00  1.8553E+00  5.2050E+00  9.2182E+00 -2.3706E+00  1.4025E+00  1.8120E+00 -3.5965E+00
             3.8778E+00

0ITERATION NO.:   52    OBJECTIVE VALUE:  -3919.52279954998        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1510
 NPARAMETR:  9.9672E-01  1.0398E+00  1.0512E+00  1.0065E+00  9.6577E-01  1.0385E+00  9.2881E-01  8.8277E-01  1.0664E+00  9.4185E-01
             1.0221E+00
 PARAMETER:  9.6710E-02  1.3907E-01  1.4997E-01  1.0651E-01  6.5169E-02  1.3763E-01  2.6163E-02 -2.4688E-02  1.6432E-01  4.0104E-02
             1.2181E-01
 GRADIENT:  -3.3722E+01  1.9698E+06  7.1247E+01 -4.6116E+01  3.6960E+01 -1.2622E-01  2.7393E+06 -2.2132E+02  8.3353E+05  2.7393E+06
            -2.2493E+06
 NUMSIGDIG:         8.5         3.3         8.0         8.3         8.5         2.3         3.3         7.7         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1510
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7912E-02 -2.3125E-02 -1.6488E-02  1.6879E-02 -2.6278E-02
 SE:             2.9858E-02  2.2583E-02  1.6003E-02  2.7883E-02  2.5022E-02
 N:                     100         100         100         100         100

 P VAL.:         5.4857E-01  3.0583E-01  3.0286E-01  5.4494E-01  2.9363E-01

 ETASHRINKSD(%)  1.0000E-10  2.4345E+01  4.6388E+01  6.5890E+00  1.6174E+01
 ETASHRINKVR(%)  1.0000E-10  4.2764E+01  7.1257E+01  1.2744E+01  2.9733E+01
 EBVSHRINKSD(%)  2.5599E-01  2.5395E+01  4.6719E+01  7.1647E+00  1.6968E+01
 EBVSHRINKVR(%)  5.1133E-01  4.4340E+01  7.1612E+01  1.3816E+01  3.1058E+01
 RELATIVEINF(%)  9.9487E+01  2.8126E+01  1.7980E+01  6.2563E+01  2.8117E+01
 EPSSHRINKSD(%)  2.1079E+01
 EPSSHRINKVR(%)  3.7715E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3919.5227995499836     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2265.4334397815728     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    41.34
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3919.523       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.97E-01  1.04E+00  1.05E+00  1.01E+00  9.66E-01  1.04E+00  9.29E-01  8.83E-01  1.07E+00  9.42E-01  1.02E+00
 


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
+        1.38E+10
 
 TH 2
+        6.76E+03  3.27E+09
 
 TH 3
+       -4.36E+09  1.30E+05  5.51E+09
 
 TH 4
+       -5.93E-01 -4.42E+09 -3.93E+01  1.19E+10
 
 TH 5
+       -4.20E+00  6.55E+04 -2.80E+02  6.61E+09  1.47E+10
 
 TH 6
+        3.41E+01  1.12E+04  1.02E+00 -2.49E+00 -7.91E-01  1.82E+02
 
 TH 7
+        3.78E+00 -6.24E+04 -1.37E+01  2.41E+00 -3.36E+00  1.74E+04  7.94E+09
 
 TH 8
+        7.78E+09 -4.36E+05 -4.92E+09  3.84E+01 -1.41E+01 -5.90E+00  8.35E+09  1.76E+10
 
 TH 9
+       -3.92E+09 -8.44E+04  1.07E+05  2.68E+01  5.74E+00  9.21E+03 -5.15E+04  4.43E+09  2.23E+09
 
 TH10
+       -7.30E+09 -6.39E+04  4.61E+09  5.53E+00 -7.47E+01  1.71E+04 -9.58E+04 -6.70E+05 -4.15E+09  7.72E+09
 
 TH11
+       -7.87E+03 -1.26E+06  3.49E+09 -1.98E+00 -1.11E+01 -1.30E+04  7.25E+04 -6.23E+09  3.14E+09  5.84E+09  4.42E+09
 
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
 #CPUT: Total CPU Time in Seconds,       53.840
Stop Time:
Sat Sep 18 06:21:40 CDT 2021
