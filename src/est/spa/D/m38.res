Wed Sep 29 19:59:55 CDT 2021
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
$DATA ../../../../data/spa/D/dat38.csv ignore=@
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
Current Date:       29 SEP 2021
Days until program expires : 200
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   9272.59574335858        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.0944E+02  1.7505E+02  2.7702E+01  1.0897E+02  1.5981E+02 -1.7814E+03 -7.3489E+02 -8.5209E+01 -9.9423E+02 -6.0570E+02
            -1.7434E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -661.422494058512        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3032E+00  1.1687E+00  9.0748E-01  1.7682E+00  1.3430E+00  2.0217E+00  1.5200E+00  9.3967E-01  1.3899E+00  1.3607E+00
             1.4103E+01
 PARAMETER:  3.6484E-01  2.5585E-01  2.9131E-03  6.6999E-01  3.9488E-01  8.0396E-01  5.1872E-01  3.7775E-02  4.2926E-01  4.0799E-01
             2.7464E+00
 GRADIENT:  -2.4123E+01  2.5076E+01 -1.4054E+01  5.8909E+01 -8.8413E+00  4.6113E+01  2.8199E+00  5.1559E+00  1.4872E+01  3.4153E+00
             1.8257E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -679.552308211147        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.2879E+00  1.3120E+00  2.0888E+00  1.7597E+00  4.9794E+00  1.8133E+00  2.6207E+00  4.0938E-01  1.2872E+00  1.0279E+01
             1.2247E+01
 PARAMETER:  3.5301E-01  3.7154E-01  8.3657E-01  6.6513E-01  1.7053E+00  6.9515E-01  1.0634E+00 -7.9312E-01  3.5244E-01  2.4301E+00
             2.6053E+00
 GRADIENT:   5.5412E+00  5.1479E+01  3.6951E+00  1.0048E+02 -5.6640E+00  1.7137E+01  5.0739E+00 -8.0808E-04  5.3459E-01  1.4353E+01
             9.3936E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -716.426401609129        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  1.0823E+00  8.4545E-01  9.2459E-01  1.1828E+00  4.6304E+00  1.5004E+00  1.7149E+00  7.5624E-02  7.9626E-01  7.0368E+00
             1.0322E+01
 PARAMETER:  1.7913E-01 -6.7891E-02  2.1598E-02  2.6792E-01  1.6326E+00  5.0573E-01  6.3938E-01 -2.4820E+00 -1.2783E-01  2.0512E+00
             2.4343E+00
 GRADIENT:  -9.2282E+00 -6.1589E+00  2.3562E+00 -4.7982E+01 -6.1668E+00  6.6442E+00  5.8519E+00  7.2135E-03  7.0149E+00 -3.8373E+00
             5.0390E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -725.005484735489        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:      367             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0730E+00  5.6104E-01  8.5977E-01  1.3838E+00  6.2126E+00  1.4851E+00  1.4389E+00  2.1167E-02  6.0260E-01  8.5935E+00
             9.5007E+00
 PARAMETER:  1.7045E-01 -4.7796E-01 -5.1091E-02  4.2486E-01  1.9266E+00  4.9551E-01  4.6390E-01 -3.7553E+00 -4.0651E-01  2.2510E+00
             2.3514E+00
 GRADIENT:  -1.0538E+01  8.7216E+00 -3.4372E+00  6.6823E-02  2.0009E+01  1.7920E+01  1.9003E+00 -1.1836E-04  4.4381E+00 -4.3748E+01
            -2.6419E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -727.362027501148        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      556
 NPARAMETR:  1.0728E+00  5.5685E-01  8.5798E-01  1.3818E+00  6.3250E+00  1.4986E+00  1.4421E+00  1.0036E-02  6.0368E-01  8.6104E+00
             1.0124E+01
 PARAMETER:  1.7025E-01 -4.8546E-01 -5.3178E-02  4.2337E-01  1.9445E+00  5.0452E-01  4.6609E-01 -4.5016E+00 -4.0472E-01  2.2530E+00
             2.4149E+00
 GRADIENT:  -1.9932E+01  9.7286E+00  3.9266E+00 -4.9249E+00 -3.3370E+00 -1.2798E+01  5.2793E-01  1.9938E-05  1.3208E+00 -2.4030E+00
            -3.5281E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -729.498207136505        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      735
 NPARAMETR:  1.0877E+00  4.9804E-01  7.9533E-01  1.4248E+00  8.3321E+00  1.6783E+00  1.4383E+00  1.0000E-02  5.9994E-01  9.7005E+00
             1.0453E+01
 PARAMETER:  1.8406E-01 -5.9708E-01 -1.2900E-01  4.5406E-01  2.2201E+00  6.1777E-01  4.6344E-01 -2.0661E+01 -4.1092E-01  2.3722E+00
             2.4468E+00
 GRADIENT:  -9.7331E+00  6.8314E+00  3.3355E+00 -1.5969E+01 -3.3201E+00  1.6438E+01  8.5606E-01  0.0000E+00  2.5875E+00  1.9065E+00
             2.8308E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -731.057921434824        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      912
 NPARAMETR:  1.0945E+00  4.6262E-01  7.6332E-01  1.4611E+00  9.9700E+00  1.6322E+00  1.4342E+00  1.0000E-02  5.9688E-01  9.7311E+00
             1.0018E+01
 PARAMETER:  1.9026E-01 -6.7085E-01 -1.7008E-01  4.7917E-01  2.3996E+00  5.8995E-01  4.6062E-01 -2.4222E+01 -4.1604E-01  2.3753E+00
             2.4044E+00
 GRADIENT:  -3.5189E+00  1.0990E+01  2.4383E+00  1.8306E+01  2.3696E+00  2.9505E+00  5.5865E-01  0.0000E+00  9.3199E-01 -6.7985E+00
            -1.0333E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -734.500445075535        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     1069
 NPARAMETR:  1.0686E+00  2.8996E-01  7.1235E-01  1.4517E+00  9.6453E+00  1.5743E+00  4.5449E-01  1.0000E-02  3.1898E-02  9.9335E+00
             1.0280E+01
 PARAMETER:  1.6634E-01 -1.1380E+00 -2.3918E-01  4.7274E-01  2.3665E+00  5.5378E-01 -6.8857E-01 -2.5740E+01 -3.3452E+00  2.3959E+00
             2.4302E+00
 GRADIENT:   2.3278E+01  2.5819E-01  7.7736E+00 -2.6763E+01  1.1132E+01  2.5390E+01  4.2906E-02  0.0000E+00  2.5394E-02  1.5700E+01
             1.7315E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -735.574813707325        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     1145
 NPARAMETR:  1.0561E+00  2.6273E-01  7.0766E-01  1.4603E+00  9.2668E+00  1.5359E+00  3.7112E-01  1.0000E-02  1.9130E-02  9.8131E+00
             1.0241E+01
 PARAMETER:  1.5454E-01 -1.2366E+00 -2.4579E-01  4.7863E-01  2.3264E+00  5.2909E-01 -8.9124E-01 -2.5740E+01 -3.8565E+00  2.3837E+00
             2.4264E+00
 GRADIENT:   2.0646E+01  2.4421E+01 -8.9784E+01  1.9827E+02 -3.2340E+01  1.8086E+01 -1.4686E-01  0.0000E+00 -4.9302E-02  3.3500E+01
            -3.9909E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -736.570485154062        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1309
 NPARAMETR:  1.0511E+00  2.5171E-01  7.1298E-01  1.4603E+00  9.3166E+00  1.5649E+00  7.1259E-02  1.0000E-02  1.0000E-02  9.4104E+00
             1.0292E+01
 PARAMETER:  1.4984E-01 -1.2795E+00 -2.3831E-01  4.7865E-01  2.3318E+00  5.4780E-01 -2.5414E+00 -2.5740E+01 -6.0635E+00  2.3418E+00
             2.4314E+00
 GRADIENT:   1.3996E+00  6.4961E-01  1.9158E+01 -2.9889E+01  1.5412E+01  1.0255E+01  7.9973E-04  0.0000E+00  0.0000E+00 -4.5760E+00
             4.0819E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -736.704826220944        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     1410
 NPARAMETR:  1.0502E+00  2.4965E-01  7.1385E-01  1.4611E+00  9.2674E+00  1.5706E+00  5.0792E-02  1.0000E-02  1.0000E-02  9.3795E+00
             1.0215E+01
 PARAMETER:  1.4899E-01 -1.2877E+00 -2.3708E-01  4.7916E-01  2.3265E+00  5.5148E-01 -2.8800E+00 -2.5740E+01 -6.5388E+00  2.3385E+00
             2.4238E+00
 GRADIENT:  -6.1523E+00  5.4583E+00  9.4482E+00 -1.2809E+01  4.9623E-01 -5.2999E+00  1.7304E-04  0.0000E+00  0.0000E+00 -1.1682E+01
             3.4679E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1410
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.1375E-03 -1.6638E-04 -1.7915E-04 -2.3118E-04 -4.3795E-02
 SE:             2.8008E-02  1.0631E-04  6.1311E-05  2.9805E-04  8.9212E-03
 N:                     100         100         100         100         100

 P VAL.:         7.4424E-01  1.1757E-01  3.4784E-03  4.3797E-01  9.1619E-07

 ETASHRINKSD(%)  6.1688E+00  9.9644E+01  9.9795E+01  9.9001E+01  7.0113E+01
 ETASHRINKVR(%)  1.1957E+01  9.9999E+01  1.0000E+02  9.9990E+01  9.1067E+01
 EBVSHRINKSD(%)  4.7289E+00  9.9665E+01  9.9761E+01  9.9054E+01  6.4921E+01
 EBVSHRINKVR(%)  9.2342E+00  9.9999E+01  9.9999E+01  9.9991E+01  8.7694E+01
 RELATIVEINF(%)  4.8156E+01  3.4597E-05  1.0691E-04  1.8562E-04  6.0161E+00
 EPSSHRINKSD(%)  5.5138E+00
 EPSSHRINKVR(%)  1.0724E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -736.70482622094357     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1.5539996572053951     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.33
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.98
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -736.705       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  2.50E-01  7.14E-01  1.46E+00  9.27E+00  1.57E+00  5.08E-02  1.00E-02  1.00E-02  9.38E+00  1.02E+01
 


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
+        2.38E+04
 
 TH 2
+        1.70E+00  5.98E+03
 
 TH 3
+       -1.07E+02 -5.37E+02  2.07E+04
 
 TH 4
+       -1.23E+02  2.56E+02 -1.62E+02  1.76E+03
 
 TH 5
+       -2.84E+00 -1.99E+02  1.24E+01 -1.04E+02  4.56E+00
 
 TH 6
+       -4.22E+01 -5.20E+01  8.94E+01 -8.46E+01  3.69E+00  8.00E+02
 
 TH 7
+        5.22E-04  5.93E-04 -1.55E-01  1.67E-02 -8.13E-04 -4.16E-02 -1.08E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -1.18E+00  2.01E+00 -9.34E+00  9.65E-01 -4.66E+00  1.67E+01  3.87E-03  0.00E+00  0.00E+00  8.70E+00
 
 TH11
+       -2.12E+01 -3.95E+00  5.59E+00 -7.51E+01  4.94E+00 -1.20E+01  5.00E-04  0.00E+00  0.00E+00 -6.37E+00  1.18E+01
 
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
 #CPUT: Total CPU Time in Seconds,       37.387
Stop Time:
Wed Sep 29 20:00:34 CDT 2021
