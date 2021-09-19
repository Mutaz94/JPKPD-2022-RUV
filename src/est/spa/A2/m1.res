Sat Sep 18 09:35:06 CDT 2021
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
$DATA ../../../../data/spa/A2/dat1.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m1.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -984.363654087131        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0937E+01  1.8236E+01  7.6193E+01 -6.7449E+01  3.3686E+01  9.6675E+00 -4.3056E+01 -3.0005E+01 -7.4874E+01 -5.7200E+01
            -1.1503E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1379.07483413681        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0355E+00  9.7697E-01  9.4614E-01  1.1102E+00  9.5201E-01  8.8757E-01  1.0990E+00  1.0285E+00  1.1586E+00  9.0847E-01
             2.9446E+00
 PARAMETER:  1.3492E-01  7.6698E-02  4.4638E-02  2.0452E-01  5.0819E-02 -1.9262E-02  1.9441E-01  1.2810E-01  2.4721E-01  4.0034E-03
             1.1800E+00
 GRADIENT:   4.9391E+01  1.4479E+01 -2.1299E+00  2.8184E+01 -1.1450E+01 -2.2164E+01  3.7021E+00  2.7634E+00  1.4577E+01  1.5003E+01
             6.6394E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1383.67173010261        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0357E+00  7.9593E-01  6.0909E-01  1.2001E+00  6.6450E-01  9.3752E-01  1.5505E+00  7.0637E-01  1.0258E+00  5.5418E-01
             2.7858E+00
 PARAMETER:  1.3510E-01 -1.2825E-01 -3.9578E-01  2.8242E-01 -3.0871E-01  3.5479E-02  5.3859E-01 -2.4762E-01  1.2551E-01 -4.9026E-01
             1.1245E+00
 GRADIENT:   4.2180E+01  2.6491E+01 -2.1507E+01  6.3600E+01  8.6737E+00 -6.2470E+00  1.4752E+01  4.2506E+00  8.9748E+00  6.9859E+00
             5.4205E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1392.09978858537        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0040E+00  7.7109E-01  2.9981E-01  1.0399E+00  4.2461E-01  9.7521E-01  1.3740E+00  3.3863E-01  1.0036E+00  3.1551E-01
             2.2891E+00
 PARAMETER:  1.0395E-01 -1.5995E-01 -1.1046E+00  1.3909E-01 -7.5658E-01  7.4897E-02  4.1776E-01 -9.8285E-01  1.0357E-01 -1.0536E+00
             9.2815E-01
 GRADIENT:  -1.0621E+01  6.6148E+01  6.7836E+01 -8.9499E+00 -1.0085E+02  1.9269E+00  4.9283E+00 -1.6796E+00  3.8232E+00 -1.8942E+00
            -1.2194E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1398.32587399707        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.8749E-01  6.1538E-01  2.0070E-01  1.0006E+00  3.3157E-01  9.7614E-01  1.4096E+00  2.5154E-01  9.7433E-01  2.2253E-01
             2.1174E+00
 PARAMETER:  8.7414E-02 -3.8552E-01 -1.5059E+00  1.0059E-01 -1.0039E+00  7.5847E-02  4.4329E-01 -1.2802E+00  7.3992E-02 -1.4027E+00
             8.5021E-01
 GRADIENT:  -2.6649E+01 -2.6935E+00  1.1541E+01 -5.3092E+00  5.0264E+00 -3.3813E+00  6.3197E+00 -2.4024E+00 -1.2540E+01 -3.1551E+00
            -2.0804E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1403.12918830192        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0018E+00  5.7929E-01  1.8001E-01  1.0126E+00  3.0032E-01  9.9585E-01  1.2141E+00  1.1454E+00  1.0618E+00  4.3458E-01
             2.0657E+00
 PARAMETER:  1.0184E-01 -4.4596E-01 -1.6147E+00  1.1256E-01 -1.1029E+00  9.5846E-02  2.9399E-01  2.3575E-01  1.6000E-01 -7.3337E-01
             8.2548E-01
 GRADIENT:   1.0358E+01  4.6908E+00  1.8122E+00  2.7594E+01  1.6705E+01  7.1963E+00  4.6921E+00  5.1516E+00 -3.8570E+00  9.2095E+00
             3.2321E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1406.50590769692        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      470
 NPARAMETR:  9.8472E-01  4.8055E-01  1.4576E-01  9.8718E-01  2.5134E-01  9.8457E-01  1.2762E+00  1.2585E+00  1.1348E+00  3.4382E-01
             1.7637E+00
 PARAMETER:  8.4605E-02 -6.3283E-01 -1.8258E+00  8.7095E-02 -1.2809E+00  8.4449E-02  3.4390E-01  3.2989E-01  2.2647E-01 -9.6763E-01
             6.6739E-01
 GRADIENT:  -2.4320E+01  1.2410E+01 -1.4523E+01  1.8271E+01 -2.2178E+01 -2.5211E-02  1.9045E+00 -6.4991E-02 -1.7657E+01  2.1629E-01
             4.0934E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1408.38759761882        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      646
 NPARAMETR:  9.9562E-01  5.5897E-01  1.6900E-01  9.7721E-01  2.8708E-01  9.8042E-01  1.2394E+00  1.3269E+00  1.1511E+00  2.3912E-01
             1.7862E+00
 PARAMETER:  9.5613E-02 -4.8166E-01 -1.6778E+00  7.6945E-02 -1.1480E+00  8.0222E-02  3.1460E-01  3.8287E-01  2.4073E-01 -1.3308E+00
             6.8010E-01
 GRADIENT:   2.3241E+00  5.3318E+00  1.0641E+01 -1.1201E+01 -1.4674E+01  1.4646E+00  4.4865E-01  4.1521E-01  2.0520E+00  1.0264E+00
             2.2931E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1408.85913036512        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      821
 NPARAMETR:  9.9300E-01  6.1149E-01  1.7204E-01  9.7185E-01  3.0702E-01  9.7322E-01  1.2187E+00  1.3366E+00  1.1241E+00  1.0406E-01
             1.8067E+00
 PARAMETER:  9.2978E-02 -3.9186E-01 -1.6600E+00  7.1446E-02 -1.0808E+00  7.2857E-02  2.9777E-01  3.9012E-01  2.1697E-01 -2.1628E+00
             6.9151E-01
 GRADIENT:   2.5973E-01 -5.6806E-01 -4.2108E-01  2.7028E-01  2.0913E+00 -2.2164E-01 -2.2574E-01 -3.5863E-01  3.7420E-01  1.3799E-01
            -2.6860E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1408.91900814143        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      996
 NPARAMETR:  9.9306E-01  6.0844E-01  1.7196E-01  9.7236E-01  3.0600E-01  9.7391E-01  1.2307E+00  1.3448E+00  1.1249E+00  2.1243E-02
             1.8066E+00
 PARAMETER:  9.3036E-02 -3.9685E-01 -1.6605E+00  7.1968E-02 -1.0842E+00  7.3561E-02  3.0760E-01  3.9622E-01  2.1769E-01 -3.7517E+00
             6.9142E-01
 GRADIENT:   1.0897E-01  1.6361E-01  4.4333E-02 -7.6849E-04 -2.8400E-02  2.9645E-02  2.2837E-01  2.5675E-02  2.1036E-01  4.0939E-03
             2.0433E-02

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1408.92074552799        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1123
 NPARAMETR:  9.9300E-01  6.0865E-01  1.7194E-01  9.7222E-01  3.0607E-01  9.7381E-01  1.2292E+00  1.3444E+00  1.1241E+00  1.0000E-02
             1.8069E+00
 PARAMETER:  9.2978E-02 -3.9652E-01 -1.6606E+00  7.1827E-02 -1.0840E+00  7.3463E-02  3.0640E-01  3.9592E-01  2.1697E-01 -4.6798E+00
             6.9163E-01
 GRADIENT:  -6.5132E-03  2.7207E-03 -1.1607E-02 -1.2834E-02  1.0386E-03 -3.7551E-03 -1.3131E-02 -3.3778E-03 -1.1691E-02  0.0000E+00
            -6.9289E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1123
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -7.5186E-04  9.8329E-03 -1.8359E-02 -1.0581E-02  2.3561E-04
 SE:             2.9462E-02  2.4275E-02  1.9666E-02  2.6819E-02  3.2922E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7964E-01  6.8544E-01  3.5054E-01  6.9319E-01  4.7420E-01

 ETASHRINKSD(%)  1.2989E+00  1.8674E+01  3.4116E+01  1.0152E+01  9.8897E+01
 ETASHRINKVR(%)  2.5809E+00  3.3862E+01  5.6592E+01  1.9273E+01  9.9988E+01
 EBVSHRINKSD(%)  1.5080E+00  1.8660E+01  3.4161E+01  1.0128E+01  9.8938E+01
 EBVSHRINKVR(%)  2.9933E+00  3.3837E+01  5.6652E+01  1.9229E+01  9.9989E+01
 RELATIVEINF(%)  9.4968E+01  5.4360E+00  1.5034E+01  5.0673E+01  7.0789E-04
 EPSSHRINKSD(%)  4.5357E+01
 EPSSHRINKVR(%)  7.0141E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1408.9207455279886     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -673.76991896425045     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.79
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1408.921       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.93E-01  6.09E-01  1.72E-01  9.72E-01  3.06E-01  9.74E-01  1.23E+00  1.34E+00  1.12E+00  1.00E-02  1.81E+00
 


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
+        1.16E+03
 
 TH 2
+        5.45E+00  1.35E+03
 
 TH 3
+       -2.15E+02  1.65E+03  9.39E+03
 
 TH 4
+       -2.86E+01  1.76E+02 -8.18E+02  7.58E+02
 
 TH 5
+        1.21E+02 -3.93E+03 -7.59E+03  1.39E+02  1.36E+04
 
 TH 6
+        3.28E-01 -5.51E+00  5.89E+00 -9.07E+00  3.67E+01  1.99E+02
 
 TH 7
+       -1.04E+00  4.99E+01 -7.30E+01 -1.27E+01  2.66E+01 -1.42E-01  6.00E+01
 
 TH 8
+        1.44E+00  6.33E+00 -7.44E+01 -4.14E+00 -2.79E+01  4.35E+00  3.42E+00  3.13E+01
 
 TH 9
+        5.61E+00 -1.91E+01  4.30E+01 -1.88E+00  1.84E+02 -5.86E-01  1.61E+01 -8.19E+00  1.06E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.52E+01 -1.92E+01 -3.55E+01 -1.02E+00 -5.06E+00  9.28E-01  5.72E+00  1.25E+01  8.55E+00  0.00E+00  5.39E+01
 
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
 #CPUT: Total CPU Time in Seconds,       19.410
Stop Time:
Sat Sep 18 09:35:27 CDT 2021
