Thu Sep 30 09:52:33 CDT 2021
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
$DATA ../../../../data/spa2/D/dat84.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m84.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   29210.9145754181        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.9884E+02  5.5979E+02 -3.5275E+01  5.9705E+02  6.4929E+01 -2.4841E+03 -1.3906E+03 -7.3948E+00 -1.8913E+03 -5.5794E+02
            -5.6267E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -490.940182203305        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.3172E+00  1.2861E+00  1.0341E+00  1.3592E+00  1.0792E+00  1.7525E+00  1.6271E+00  9.6740E-01  1.4166E+00  9.7207E-01
             1.4405E+01
 PARAMETER:  3.7549E-01  3.5158E-01  1.3352E-01  4.0686E-01  1.7621E-01  6.6106E-01  5.8678E-01  6.6853E-02  4.4826E-01  7.1671E-02
             2.7676E+00
 GRADIENT:   3.5423E+00 -7.4233E+00 -1.5811E+01  3.5808E+01  1.0302E+01  2.8265E+01 -2.6094E+01  3.2231E+00 -8.2093E+00  1.0344E+01
             1.1281E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -558.474698942778        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.2575E+00  9.2637E-01  5.8944E+00  1.7737E+00  1.7863E+00  2.2685E+00  5.5988E+00  5.6992E-01  2.0587E+00  4.4442E-01
             1.3157E+01
 PARAMETER:  3.2916E-01  2.3519E-02  1.8740E+00  6.7309E-01  6.8015E-01  9.1912E-01  1.8225E+00 -4.6226E-01  8.2209E-01 -7.1099E-01
             2.6770E+00
 GRADIENT:   5.9607E-01  7.3148E+00 -3.8128E+00  1.4008E+01 -3.3775E+00  5.6556E+01  4.1266E+01  1.2176E-01  2.7822E+01  1.7074E+00
             1.6174E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -564.585327057601        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  1.2454E+00  1.0614E+00  5.8813E+00  1.5615E+00  1.9677E+00  1.9229E+00  4.7202E+00  4.8093E-01  1.6360E+00  9.1168E-01
             1.3159E+01
 PARAMETER:  3.1946E-01  1.5958E-01  1.8718E+00  5.4568E-01  7.7685E-01  7.5381E-01  1.6519E+00 -6.3202E-01  5.9228E-01  7.5372E-03
             2.6771E+00
 GRADIENT:   1.8200E+00 -1.1438E+00 -2.6718E+00  7.1074E+00  3.6485E-02  2.8500E+01  7.0995E+00  6.2159E-02  1.3487E+01  6.4508E+00
             1.5330E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -571.724456874887        NO. OF FUNC. EVALS.: 154
 CUMULATIVE NO. OF FUNC. EVALS.:      386
 NPARAMETR:  1.2458E+00  1.0604E+00  5.8723E+00  1.5593E+00  1.9651E+00  1.9236E+00  4.7414E+00  1.2234E-01  1.2890E+00  2.5861E-01
             1.2793E+01
 PARAMETER:  3.1981E-01  1.5861E-01  1.8702E+00  5.4422E-01  7.7556E-01  7.5417E-01  1.6563E+00 -2.0010E+00  3.5389E-01 -1.2524E+00
             2.6489E+00
 GRADIENT:   8.2730E+00  3.0045E+00 -2.0376E+00  3.2176E+01  1.2109E+00  2.9137E+01  2.6090E-03  4.4158E-03 -2.3920E+00  5.1706E-01
             1.1573E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -576.140072923271        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      501
 NPARAMETR:  1.2455E+00  1.0602E+00  5.8823E+00  1.5591E+00  1.9661E+00  1.9208E+00  4.7364E+00  1.0000E-02  1.4779E+00  2.2960E-02
             1.1705E+01
 PARAMETER:  3.1950E-01  1.5849E-01  1.8719E+00  5.4411E-01  7.7605E-01  7.5276E-01  1.6553E+00 -4.5978E+00  4.9061E-01 -3.6740E+00
             2.5600E+00
 GRADIENT:   1.4071E+01  6.4524E+00 -2.8222E+00  3.0246E+01  4.6132E-01  1.2560E+01 -2.5242E+01  0.0000E+00 -7.3399E-01  4.3538E-03
            -2.1845E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -582.920103745081        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:      615
 NPARAMETR:  1.1881E+00  9.4956E-01  8.3319E+00  1.4084E+00  1.9566E+00  1.7821E+00  5.1636E+00  1.0000E-02  1.4035E+00  1.2438E-02
             1.1540E+01
 PARAMETER:  2.7235E-01  4.8240E-02  2.2201E+00  4.4244E-01  7.7121E-01  6.7778E-01  1.7416E+00 -4.5387E+00  4.3896E-01 -4.2870E+00
             2.5458E+00
 GRADIENT:   1.8322E+01 -1.3401E+00 -3.5020E-01 -3.2818E+00 -8.5809E-01  1.0929E+01  3.3309E+01  2.8840E-06  8.6931E+00  1.5079E-03
             3.8022E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -584.224178941792        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      796
 NPARAMETR:  1.1765E+00  8.8636E-01  1.0387E+01  1.4423E+00  2.0460E+00  1.7596E+00  5.5646E+00  1.0000E-02  1.1647E+00  1.0000E-02
             1.1533E+01
 PARAMETER:  2.6256E-01 -2.0637E-02  2.4405E+00  4.6624E-01  8.1589E-01  6.6508E-01  1.8164E+00 -4.5457E+00  2.5243E-01 -4.6082E+00
             2.5452E+00
 GRADIENT:   2.5698E-01  1.4403E-01 -2.7768E-01 -1.8009E+00 -1.8941E-01  2.1569E-01  4.6860E-01  3.9147E-07 -6.0191E-01  0.0000E+00
             1.5572E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -585.271845104097        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      974
 NPARAMETR:  1.1780E+00  5.2553E-01  4.3455E+01  1.7827E+00  2.1947E+00  1.7534E+00  6.4322E+00  1.0000E-02  1.5093E+00  1.0000E-02
             1.1576E+01
 PARAMETER:  2.6384E-01 -5.4335E-01  3.8717E+00  6.7814E-01  8.8604E-01  6.6156E-01  1.9613E+00 -4.5665E+00  5.1167E-01 -4.6790E+00
             2.5489E+00
 GRADIENT:  -1.7158E+00 -3.8811E-02 -6.4348E-03  4.4745E-01 -2.6878E-01  7.0527E-01 -1.3266E-01  0.0000E+00 -3.5888E-01  0.0000E+00
             2.0130E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -585.365172402208        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1157             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1770E+00  5.0418E-01  5.3183E+01  1.7855E+00  2.2048E+00  1.7549E+00  6.6787E+00  1.0305E-02  1.5223E+00  1.0000E-02
             1.1540E+01
 PARAMETER:  2.6294E-01 -5.8483E-01  4.0737E+00  6.7969E-01  8.9063E-01  6.6241E-01  1.9989E+00 -4.4751E+00  5.2024E-01 -4.6790E+00
             2.5458E+00
 GRADIENT:   1.0982E+01  2.2666E+00  1.7235E-03  1.2481E+01  1.0216E+00  8.8671E+00  6.3085E+01  4.4976E-07  3.3658E+00  0.0000E+00
             3.1590E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -585.421651076892        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1338             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1815E+00  4.8441E-01  5.7794E+01  1.8156E+00  2.2115E+00  1.7496E+00  6.7910E+00  1.0748E-02  1.5424E+00  1.0000E-02
             1.1561E+01
 PARAMETER:  2.6681E-01 -6.2482E-01  4.1569E+00  6.9643E-01  8.9366E-01  6.5937E-01  2.0156E+00 -4.4330E+00  5.3332E-01 -4.6790E+00
             2.5476E+00
 GRADIENT:   1.2752E+01  2.8044E+00 -3.9370E-03  1.4499E+01  8.7427E-01  7.9999E+00  6.5063E+01  4.0942E-07  3.4569E+00  0.0000E+00
             3.1626E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -585.460082808612        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1513
 NPARAMETR:  1.1835E+00  4.5990E-01  6.8714E+01  1.8460E+00  2.2209E+00  1.7533E+00  6.7352E+00  1.0734E-02  1.5528E+00  1.0000E-02
             1.1585E+01
 PARAMETER:  2.6849E-01 -6.7676E-01  4.3300E+00  7.1301E-01  8.9792E-01  6.6149E-01  2.0074E+00 -4.4344E+00  5.4009E-01 -4.6790E+00
             2.5497E+00
 GRADIENT:   5.6108E-01  2.2311E-01 -3.3985E-03  1.5861E-01 -1.4560E-01  7.6033E-01  2.4032E+00  2.8814E-07 -4.2173E-01  0.0000E+00
             1.6241E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -585.519525786844        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1693
 NPARAMETR:  1.1836E+00  4.3094E-01  8.5835E+01  1.8713E+00  2.2262E+00  1.7531E+00  6.8442E+00  1.0000E-02  1.5708E+00  1.0000E-02
             1.1584E+01
 PARAMETER:  2.6855E-01 -7.4178E-01  4.5524E+00  7.2665E-01  9.0031E-01  6.6140E-01  2.0234E+00 -4.5233E+00  5.5161E-01 -4.6790E+00
             2.5497E+00
 GRADIENT:   5.2158E-01  1.0457E-01 -2.6599E-04  6.4100E-02 -2.3173E-01  8.0563E-01  2.6034E+00  5.9315E-08 -5.4321E-01  0.0000E+00
             1.5181E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -585.577675621971        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1873
 NPARAMETR:  1.1839E+00  4.0805E-01  1.0030E+02  1.8938E+00  2.2313E+00  1.7526E+00  6.9638E+00  1.0000E-02  1.5895E+00  1.0000E-02
             1.1584E+01
 PARAMETER:  2.6886E-01 -7.9637E-01  4.7082E+00  7.3857E-01  9.0259E-01  6.6110E-01  2.0407E+00 -4.5233E+00  5.6344E-01 -4.6790E+00
             2.5497E+00
 GRADIENT:   6.0293E-01  2.1803E-01 -1.6588E-03 -1.2078E-01 -1.8434E-01  6.8215E-01  3.5283E+00  4.2101E-08 -3.6732E-01  0.0000E+00
             1.4175E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -585.650038980364        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2059             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1820E+00  3.7739E-01  1.3800E+02  1.9111E+00  2.2389E+00  1.7468E+00  7.1864E+00  2.2108E-02  1.6059E+00  1.0000E-02
             1.1548E+01
 PARAMETER:  2.6725E-01 -8.7448E-01  5.0272E+00  7.4766E-01  9.0598E-01  6.5778E-01  2.0722E+00 -3.7118E+00  5.7369E-01 -4.6790E+00
             2.5465E+00
 GRADIENT:   1.2875E+01  2.6250E+00 -1.7089E-03  1.7357E+01  9.6388E-01  7.8537E+00  7.0915E+01  4.0652E-07  3.2629E+00  0.0000E+00
             2.9605E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -585.684574608838        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     2218
 NPARAMETR:  1.1826E+00  3.6183E-01  1.9789E+02  1.9290E+00  2.2404E+00  1.7476E+00  7.2569E+00  2.3358E-02  1.6203E+00  1.0000E-02
             1.1556E+01
 PARAMETER:  2.6773E-01 -9.1658E-01  5.3877E+00  7.5702E-01  9.0667E-01  6.5823E-01  2.0820E+00 -3.6568E+00  5.8260E-01 -4.6790E+00
             2.5472E+00
 GRADIENT:   1.2855E+01  2.6796E+00  1.3670E-03  1.8285E+01  6.8696E-01  8.0113E+00  7.1744E+01  1.8606E-07  3.4592E+00  0.0000E+00
             2.9928E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -585.709166008443        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2405             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1824E+00  3.5004E-01  1.8395E+02  1.9426E+00  2.2444E+00  1.7493E+00  7.3217E+00  2.2517E-02  1.6301E+00  1.0000E-02
             1.1563E+01
 PARAMETER:  2.6757E-01 -9.4970E-01  5.3146E+00  7.6403E-01  9.0843E-01  6.5924E-01  2.0908E+00 -3.6935E+00  5.8865E-01 -4.6790E+00
             2.5479E+00
 GRADIENT:   1.2412E+01  2.7260E+00 -1.4637E-03  1.8813E+01  9.3087E-01  8.3167E+00  7.2729E+01  2.0583E-07  3.6376E+00  0.0000E+00
             3.0564E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -585.725269489907        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2581
 NPARAMETR:  1.1835E+00  3.3843E-01  2.4079E+02  1.9550E+00  2.2464E+00  1.7506E+00  7.3361E+00  2.2773E-02  1.6342E+00  1.0000E-02
             1.1574E+01
 PARAMETER:  2.6846E-01 -9.8345E-01  5.5839E+00  7.7039E-01  9.0935E-01  6.5998E-01  2.0928E+00 -3.6822E+00  5.9117E-01 -4.6790E+00
             2.5487E+00
 GRADIENT:   3.1423E-01  1.3400E-01  1.9440E-04 -7.7540E-01 -6.0156E-02  2.9785E-01  5.3697E+00  1.0313E-07 -1.3536E-01  0.0000E+00
             2.5974E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -585.760529870281        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2766             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1834E+00  3.0972E-01  3.2040E+02  1.9705E+00  2.2506E+00  1.7475E+00  7.4928E+00  2.3341E-02  1.6462E+00  1.0000E-02
             1.1562E+01
 PARAMETER:  2.6835E-01 -1.0721E+00  5.8696E+00  7.7829E-01  9.1120E-01  6.5821E-01  2.1139E+00 -3.6576E+00  5.9849E-01 -4.6790E+00
             2.5477E+00
 GRADIENT:   1.2872E+01  2.1983E+00 -4.7252E-04  1.8586E+01  1.1817E+00  8.1645E+00  7.4992E+01  7.8534E-08  3.7262E+00  0.0000E+00
             3.0145E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -585.773526839465        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2951             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1830E+00  3.0254E-01  3.1670E+02  1.9803E+00  2.2505E+00  1.7500E+00  7.5378E+00  2.2852E-02  1.6505E+00  1.0000E-02
             1.1570E+01
 PARAMETER:  2.6808E-01 -1.0955E+00  5.8580E+00  7.8325E-01  9.1117E-01  6.5961E-01  2.1199E+00 -3.6787E+00  6.0106E-01 -4.6790E+00
             2.5484E+00
 GRADIENT:   1.2381E+01  2.2386E+00 -5.0305E-04  1.9419E+01  1.1138E+00  8.6229E+00  7.5550E+01  8.7180E-08  3.6368E+00  0.0000E+00
             3.0680E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -585.787091893216        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3135             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1839E+00  2.9460E-01  3.3733E+02  1.9904E+00  2.2518E+00  1.7497E+00  7.5999E+00  2.4270E-02  1.6573E+00  1.0000E-02
             1.1575E+01
 PARAMETER:  2.6881E-01 -1.1221E+00  5.9211E+00  7.8833E-01  9.1173E-01  6.5946E-01  2.1281E+00 -3.6185E+00  6.0521E-01 -4.6790E+00
             2.5489E+00
 GRADIENT:   1.2623E+01  2.3315E+00 -6.5574E-04  1.9945E+01  1.1056E+00  8.5518E+00  7.6633E+01  7.7791E-08  3.7896E+00  0.0000E+00
             3.0772E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -585.791972713899        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     3313
 NPARAMETR:  1.1841E+00  2.9703E-01  4.4856E+02  1.9960E+00  2.2510E+00  1.7502E+00  7.6328E+00  2.2237E-02  1.6601E+00  1.0000E-02
             1.1579E+01
 PARAMETER:  2.6900E-01 -1.1139E+00  6.2060E+00  7.9112E-01  9.1135E-01  6.5975E-01  2.1325E+00 -3.7060E+00  6.0689E-01 -4.6790E+00
             2.5492E+00
 GRADIENT:   2.5542E-01  2.5921E-01  2.7173E-04 -5.0048E-01 -1.1345E-01  1.9171E-01  7.2580E+00  2.7787E-08 -4.9609E-03  0.0000E+00
             5.5614E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -585.799298743077        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     3497
 NPARAMETR:  1.1838E+00  2.9065E-01  4.1618E+02  2.0007E+00  2.2500E+00  1.7496E+00  7.6825E+00  2.1299E-02  1.6623E+00  1.0000E-02
             1.1574E+01
 PARAMETER:  2.6871E-01 -1.1356E+00  6.1311E+00  7.9350E-01  9.1095E-01  6.5939E-01  2.1389E+00 -3.7491E+00  6.0819E-01 -4.6790E+00
             2.5488E+00
 GRADIENT:   1.9731E-01  2.6836E-01  8.2573E-05 -5.7749E-01 -1.1044E-01  7.8162E-02  7.6267E+00  5.0656E-08 -7.1557E-02  0.0000E+00
            -4.6536E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -585.805535820614        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     3653
 NPARAMETR:  1.1836E+00  2.8445E-01  4.6515E+02  2.0070E+00  2.2486E+00  1.7494E+00  7.7252E+00  2.6404E-02  1.6650E+00  1.0000E-02
             1.1571E+01
 PARAMETER:  2.6856E-01 -1.1572E+00  6.2424E+00  7.9663E-01  9.1031E-01  6.5926E-01  2.1445E+00 -3.5342E+00  6.0980E-01 -4.6790E+00
             2.5485E+00
 GRADIENT:   1.3766E-01  2.7598E-01  4.8527E-04 -2.1343E-01 -2.5293E-01  3.2225E-02  7.7518E+00  7.1974E-08 -1.8717E-01  0.0000E+00
            -9.9795E-01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -585.822225984679        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     3852
 NPARAMETR:  1.1837E+00  2.6205E-01  4.7299E+02  2.0222E+00  2.2540E+00  1.7512E+00  7.8702E+00  1.0000E-02  1.6795E+00  1.0000E-02
             1.1588E+01
 PARAMETER:  2.6863E-01 -1.2392E+00  6.2591E+00  8.0418E-01  9.1271E-01  6.6028E-01  2.1631E+00 -4.6626E+00  6.1848E-01 -4.6790E+00
             2.5499E+00
 GRADIENT:  -2.2246E-01  5.4549E-02 -4.7485E-04 -1.6628E+00  2.7408E-01  3.8782E-01  8.1350E+00  0.0000E+00  4.0281E-01  0.0000E+00
             1.2060E+00

0ITERATION NO.:  125    OBJECTIVE VALUE:  -585.825994521675        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     4039
 NPARAMETR:  1.1842E+00  2.5989E-01  3.6247E+02  2.0266E+00  2.2515E+00  1.7499E+00  7.8850E+00  1.0000E-02  1.6768E+00  1.0000E-02
             1.1579E+01
 PARAMETER:  2.6908E-01 -1.2475E+00  5.9929E+00  8.0638E-01  9.1160E-01  6.5958E-01  2.1650E+00 -4.6626E+00  6.1691E-01 -4.6790E+00
             2.5492E+00
 GRADIENT:   1.5892E-01  9.1403E-02 -7.2869E-04 -5.0859E-01  1.0438E-01  1.5859E-01  8.0066E+00  0.0000E+00 -1.5370E-02  0.0000E+00
            -3.7630E-01

0ITERATION NO.:  130    OBJECTIVE VALUE:  -585.828263422233        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     4226
 NPARAMETR:  1.1843E+00  2.5498E-01  3.1879E+02  2.0294E+00  2.2544E+00  1.7496E+00  7.9226E+00  1.0000E-02  1.6785E+00  1.0000E-02
             1.1580E+01
 PARAMETER:  2.6916E-01 -1.2666E+00  5.8645E+00  8.0774E-01  9.1291E-01  6.5939E-01  2.1697E+00 -4.6626E+00  6.1790E-01 -4.6790E+00
             2.5493E+00
 GRADIENT:   1.7951E-01  5.1394E-02 -2.0741E-03 -8.8345E-01  3.7289E-01  1.1166E-01  8.1922E+00  0.0000E+00  6.2189E-02  0.0000E+00
            -1.9543E-01

0ITERATION NO.:  135    OBJECTIVE VALUE:  -585.828731574090        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:     4369
 NPARAMETR:  1.1844E+00  2.5527E-01  7.9238E+02  2.0305E+00  2.2532E+00  1.7498E+00  7.9210E+00  1.0000E-02  1.6787E+00  1.0000E-02
             1.1582E+01
 PARAMETER:  2.6924E-01 -1.2654E+00  6.7750E+00  8.0828E-01  9.1235E-01  6.5953E-01  2.1695E+00 -4.6626E+00  6.1800E-01 -4.6790E+00
             2.5494E+00
 GRADIENT:   1.5750E-01  6.1060E-02  4.1803E-04 -4.2200E-01 -5.2161E-02  1.5810E-01  8.0597E+00  0.0000E+00  4.5347E-02  0.0000E+00
            -1.7080E-01

0ITERATION NO.:  138    OBJECTIVE VALUE:  -585.829462132063        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     4467
 NPARAMETR:  1.1843E+00  2.5455E-01  6.0889E+02  2.0308E+00  2.2541E+00  1.7497E+00  7.9387E+00  1.0000E-02  1.6785E+00  1.0000E-02
             1.1580E+01
 PARAMETER:  2.6915E-01 -1.2682E+00  6.5116E+00  8.0845E-01  9.1273E-01  6.5945E-01  2.1718E+00 -4.6626E+00  6.1791E-01 -4.6790E+00
             2.5493E+00
 GRADIENT:   1.6513E-01  1.2115E-01  6.4643E-05 -6.6236E-01  6.7986E-02  1.4030E-01  8.4429E+00  0.0000E+00 -2.4866E-03  0.0000E+00
            -2.7181E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     4467
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2025E-02  6.1309E-02  1.6810E-07 -7.8815E-02 -1.9951E-05
 SE:             2.8039E-02  1.9618E-02  7.4863E-08  1.8167E-02  3.8160E-05
 N:                     100         100         100         100         100

 P VAL.:         6.6803E-01  1.7768E-03  2.4737E-02  1.4369E-05  6.0110E-01

 ETASHRINKSD(%)  6.0658E+00  3.4279E+01  1.0000E+02  3.9138E+01  9.9872E+01
 ETASHRINKVR(%)  1.1764E+01  5.6807E+01  1.0000E+02  6.2958E+01  1.0000E+02
 EBVSHRINKSD(%)  7.8486E+00  4.5543E+01  9.9999E+01  2.5907E+01  9.9761E+01
 EBVSHRINKVR(%)  1.5081E+01  7.0344E+01  1.0000E+02  4.5102E+01  9.9999E+01
 RELATIVEINF(%)  8.1428E+01  1.7368E+01  5.0321E-10  2.7026E+01  7.8777E-05
 EPSSHRINKSD(%)  7.8994E+00
 EPSSHRINKVR(%)  1.5175E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -585.82946213206253     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       516.89677771354457     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   114.43
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    12.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -585.829       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.18E+00  2.55E-01  6.09E+02  2.03E+00  2.25E+00  1.75E+00  7.94E+00  1.00E-02  1.68E+00  1.00E-02  1.16E+01
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        4.89E+02
 
 TH 2
+       -1.30E+02  1.24E+02
 
 TH 3
+       -1.35E-04  1.40E-05  8.21E-11
 
 TH 4
+       -2.06E+02  7.47E+01  9.72E-05  1.47E+02
 
 TH 5
+        2.03E+01 -9.18E+00 -1.07E-05 -1.64E+01  1.91E+00
 
 TH 6
+       -1.56E+01 -1.38E+01  1.53E-05 -2.36E+00 -5.08E-01  4.50E+01
 
 TH 7
+        2.58E+00  8.10E+00 -6.79E-06 -3.95E+00  3.37E-01 -1.45E+00  1.30E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.73E+00 -7.27E+00 -2.58E-06 -7.15E+00  8.70E-01  2.72E+00 -2.41E-01  0.00E+00  7.04E-01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.17E+00 -8.56E-01  4.77E-08 -1.56E+00  2.34E-01  3.10E+00  2.76E-02  0.00E+00  4.01E-01  0.00E+00  7.02E-01
 
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
+        2.13E+02
 
 TH 2
+       -1.45E+00  1.03E+02
 
 TH 3
+       -1.42E-05 -2.09E-05  6.59E-10
 
 TH 4
+       -1.27E+01  1.78E+01  3.75E-05  7.14E+01
 
 TH 5
+       -2.16E+00 -4.66E+00 -9.97E-05 -7.71E+00  1.21E+01
 
 TH 6
+       -2.63E+00 -1.97E+00  2.81E-06 -3.19E+00 -4.73E-01  4.62E+01
 
 TH 7
+        9.03E-01  1.06E+01 -5.66E-06 -3.93E+00  1.13E-01  2.71E-01  1.70E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.72E+00 -4.24E+00  1.79E-06 -8.67E+00  1.77E+00 -5.11E+00 -6.06E-01  0.00E+00  2.55E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.29E+00 -1.87E+00  1.18E-06 -5.02E+00  5.26E-01  1.93E+00  1.24E-01  0.00E+00  1.86E+00  0.00E+00  4.52E+00
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.21E+02
 
 TH 2
+        6.89E+01  1.28E+02
 
 TH 3
+        2.08E-06 -5.18E-06  3.13E-10
 
 TH 4
+        8.28E+01  1.51E+01  3.69E-05  7.28E+01
 
 TH 5
+       -1.12E+01  1.32E+00 -2.72E-05 -9.75E+00  4.03E+00
 
 TH 6
+        1.71E+01  1.92E+01 -1.08E-06 -3.62E+00 -1.12E+00  4.32E+01
 
 TH 7
+       -5.45E-01  1.15E+01 -2.93E-06 -3.60E+00  1.09E+00  1.16E+00  1.77E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -7.23E+00 -5.95E+00  1.42E-05 -4.56E-01  1.58E-01  5.31E+00 -7.23E-01  0.00E+00  1.54E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.43E+01 -1.11E+01 -1.59E-05 -2.01E+01  3.82E+00 -2.90E+00  5.94E-01  0.00E+00 -1.74E+00  0.00E+00  6.64E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,      126.831
Stop Time:
Thu Sep 30 09:54:42 CDT 2021
