Wed Sep 29 18:06:43 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat35.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1703.94310234915        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4635E+02  8.8596E+00 -4.3084E+01  8.2935E+01  7.4385E+01  6.5034E+01 -8.5720E-04  1.1755E+01  6.3137E+00 -1.4018E+01
             1.2009E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1712.00610080044        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0542E+00  1.0228E+00  1.1413E+00  1.0160E+00  1.0222E+00  9.4043E-01  1.0187E+00  8.6715E-01  1.0098E+00  1.1359E+00
             9.5419E-01
 PARAMETER:  1.5276E-01  1.2252E-01  2.3213E-01  1.1590E-01  1.2194E-01  3.8577E-02  1.1852E-01 -4.2541E-02  1.0976E-01  2.2741E-01
             5.3104E-02
 GRADIENT:  -2.7895E+00  2.1756E+01 -2.5418E-01  2.4102E+01 -8.5649E+00 -1.1742E+01  1.3036E+00  3.5580E+00 -9.7929E-02 -3.7207E+00
            -8.0418E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1712.59187259670        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0537E+00  8.3681E-01  1.2627E+00  1.1134E+00  9.9302E-01  9.6431E-01  9.4244E-01  6.5294E-01  9.5620E-01  1.2406E+00
             9.7652E-01
 PARAMETER:  1.5230E-01 -7.8161E-02  3.3326E-01  2.0742E-01  9.2999E-02  6.3660E-02  4.0714E-02 -3.2626E-01  5.5212E-02  3.1560E-01
             7.6238E-02
 GRADIENT:   6.7148E-01  6.6223E+00  1.6232E+01 -9.6401E+00 -2.3662E+01 -4.9048E-01 -9.1182E-01 -3.5810E+00 -4.6850E+00  6.3552E+00
             4.3039E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1713.49554370494        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.0549E+00  7.7556E-01  1.0547E+00  1.1538E+00  8.9863E-01  9.6563E-01  1.4028E+00  5.6567E-01  8.6561E-01  1.0538E+00
             9.7264E-01
 PARAMETER:  1.5342E-01 -1.5416E-01  1.5323E-01  2.4303E-01 -6.8885E-03  6.5028E-02  4.3849E-01 -4.6975E-01 -4.4324E-02  1.5245E-01
             7.2259E-02
 GRADIENT:   8.9446E-01  7.5891E+00 -4.3787E+00  1.5018E+01  4.2426E+00 -5.7397E-01  1.0121E+00  8.7940E-01  1.5997E+00 -2.1500E-02
             4.4413E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1714.03799920969        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      724
 NPARAMETR:  1.0531E+00  5.9620E-01  1.0078E+00  1.2475E+00  8.0818E-01  9.6396E-01  1.7235E+00  4.5710E-01  7.7493E-01  9.9337E-01
             9.7741E-01
 PARAMETER:  1.5173E-01 -4.1718E-01  1.0777E-01  3.2115E-01 -1.1297E-01  6.3290E-02  6.4434E-01 -6.8286E-01 -1.5499E-01  9.3351E-02
             7.7155E-02
 GRADIENT:  -1.7159E-01  4.4562E+00  3.9225E+00  2.9953E+00 -8.4684E+00 -5.9584E-01 -2.2613E+00  2.0074E-02 -4.5255E+00 -4.3669E-01
             1.4511E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1714.24812175448        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      901
 NPARAMETR:  1.0503E+00  4.3439E-01  1.0264E+00  1.3405E+00  7.7194E-01  9.6304E-01  2.1891E+00  3.7967E-01  7.4660E-01  1.0052E+00
             9.7954E-01
 PARAMETER:  1.4911E-01 -7.3381E-01  1.2603E-01  3.9307E-01 -1.5885E-01  6.2335E-02  8.8351E-01 -8.6845E-01 -1.9223E-01  1.0522E-01
             7.9328E-02
 GRADIENT:  -7.1745E-01  3.0709E+00  4.4550E+00  3.3562E+00 -6.3465E+00 -2.0049E-02  9.9346E-02 -9.4841E-01 -1.6413E-01 -2.9799E-01
            -1.1497E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1714.43033386270        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1084
 NPARAMETR:  1.0502E+00  4.2087E-01  1.0400E+00  1.3468E+00  7.7645E-01  9.6289E-01  2.2134E+00  5.1633E-01  7.4523E-01  9.8643E-01
             9.7719E-01
 PARAMETER:  1.4898E-01 -7.6544E-01  1.3927E-01  3.9776E-01 -1.5302E-01  6.2186E-02  8.9454E-01 -5.6101E-01 -1.9406E-01  8.6335E-02
             7.6926E-02
 GRADIENT:  -2.2889E-01  4.2557E-01 -4.3432E+00 -1.3612E+00  2.0906E+00  3.9702E-02 -8.3481E-01  2.1974E-01  3.3437E-01 -1.2590E+00
             5.4113E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1714.57846261943        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1260
 NPARAMETR:  1.0469E+00  2.9803E-01  1.2210E+00  1.4431E+00  8.2397E-01  9.5995E-01  2.6916E+00  7.7337E-01  7.2369E-01  1.0681E+00
             9.8557E-01
 PARAMETER:  1.4586E-01 -1.1106E+00  2.9969E-01  4.6677E-01 -9.3622E-02  5.9121E-02  1.0901E+00 -1.5700E-01 -2.2339E-01  1.6588E-01
             8.5464E-02
 GRADIENT:  -6.9992E-01  2.2991E+00 -1.2359E+01  1.5748E+01  6.5584E+00 -4.6131E-02 -1.0456E+00  2.4039E+00  1.4035E+00  2.0082E+00
             5.6312E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1714.58375476416        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1435
 NPARAMETR:  1.0456E+00  2.4399E-01  1.3296E+00  1.4883E+00  8.5259E-01  9.5870E-01  2.9920E+00  8.8143E-01  7.1461E-01  1.1154E+00
             9.8727E-01
 PARAMETER:  1.4461E-01 -1.3106E+00  3.8491E-01  4.9762E-01 -5.9471E-02  5.7821E-02  1.1960E+00 -2.6214E-02 -2.3602E-01  2.0923E-01
             8.7193E-02
 GRADIENT:  -8.6931E-01  3.2286E+00 -1.2721E+01  2.5273E+01  6.2116E+00 -1.5123E-01 -6.0268E-01  2.6718E+00  1.0089E+00  3.2324E+00
             6.8031E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1714.59911131759        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1610
 NPARAMETR:  1.0448E+00  2.0463E-01  1.4333E+00  1.5217E+00  8.7947E-01  9.5791E-01  3.2671E+00  9.7067E-01  7.0834E-01  1.1549E+00
             9.8696E-01
 PARAMETER:  1.4379E-01 -1.4866E+00  4.5997E-01  5.1982E-01 -2.8439E-02  5.7003E-02  1.2839E+00  7.0226E-02 -2.4483E-01  2.4400E-01
             8.6876E-02
 GRADIENT:  -9.3279E-01  3.4795E+00 -1.1235E+01  2.9922E+01  5.1230E+00 -2.1797E-01 -8.2469E-02  2.3942E+00  5.5362E-01  3.8080E+00
             6.8691E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1714.64096901967        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1785
 NPARAMETR:  1.0440E+00  1.6059E-01  1.5953E+00  1.5585E+00  9.1937E-01  9.5727E-01  3.6536E+00  1.1050E+00  7.0107E-01  1.2029E+00
             9.8308E-01
 PARAMETER:  1.4302E-01 -1.7289E+00  5.6708E-01  5.4373E-01  1.5930E-02  5.6335E-02  1.3957E+00  1.9988E-01 -2.5515E-01  2.8477E-01
             8.2930E-02
 GRADIENT:  -8.6583E-01  2.9380E+00 -7.7420E+00  2.8082E+01  3.1439E+00 -2.5122E-01  5.9025E-01  1.7719E+00 -5.4270E-02  3.7016E+00
             5.3378E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1714.68750787599        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1960
 NPARAMETR:  1.0437E+00  1.3298E-01  1.7409E+00  1.5792E+00  9.5158E-01  9.5724E-01  3.9117E+00  1.2171E+00  6.9716E-01  1.2307E+00
             9.7714E-01
 PARAMETER:  1.4280E-01 -1.9175E+00  6.5439E-01  5.5689E-01  5.0370E-02  5.6294E-02  1.4640E+00  2.9647E-01 -2.6074E-01  3.0757E-01
             7.6875E-02
 GRADIENT:  -3.5750E-01  1.7242E+00 -3.6729E+00  1.4826E+01  1.6467E+00 -1.5158E-01  1.0948E+00  8.0247E-01 -4.8726E-01  2.0722E+00
             2.6562E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1714.87595901182        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2151
 NPARAMETR:  1.0454E+00  1.2743E-01  1.7511E+00  1.5702E+00  9.4909E-01  9.5785E-01  3.7918E+00  1.2264E+00  7.0065E-01  1.2221E+00
             9.7181E-01
 PARAMETER:  1.4437E-01 -1.9602E+00  6.6024E-01  5.5117E-01  4.7744E-02  5.6934E-02  1.4329E+00  3.0406E-01 -2.5575E-01  3.0054E-01
             7.1409E-02
 GRADIENT:   4.0722E+00  6.5171E-01 -8.1186E-01 -2.1623E+01  4.9645E-01  1.2860E-01  1.8278E+00  5.0701E-01 -1.2786E+00  2.2486E-01
             2.4608E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1714.88197087191        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2334
 NPARAMETR:  1.0453E+00  1.2650E-01  1.7524E+00  1.5705E+00  9.4895E-01  9.5787E-01  3.7884E+00  1.2192E+00  7.0260E-01  1.2212E+00
             9.7147E-01
 PARAMETER:  1.4429E-01 -1.9675E+00  6.6098E-01  5.5137E-01  4.7600E-02  5.6953E-02  1.4319E+00  2.9821E-01 -2.5296E-01  2.9984E-01
             7.1059E-02
 GRADIENT:   3.9149E+00  4.6126E-01 -4.6020E-02 -2.1460E+01  2.2769E-01  1.3818E-01  1.1139E+00  5.8238E-02 -4.9715E-02  1.5760E-02
             5.0015E-02

0ITERATION NO.:   66    OBJECTIVE VALUE:  -1714.88197087191        NO. OF FUNC. EVALS.:  24
 CUMULATIVE NO. OF FUNC. EVALS.:     2358
 NPARAMETR:  1.0453E+00  1.2650E-01  1.7524E+00  1.5705E+00  9.4895E-01  9.5787E-01  3.7884E+00  1.2192E+00  7.0260E-01  1.2212E+00
             9.7147E-01
 PARAMETER:  1.4429E-01 -1.9675E+00  6.6098E-01  5.5137E-01  4.7600E-02  5.6953E-02  1.4319E+00  2.9821E-01 -2.5296E-01  2.9984E-01
             7.1059E-02
 GRADIENT:   2.3632E-02 -1.0141E-02 -2.0000E-03 -1.1470E-01  2.3694E-01 -6.5498E-04 -5.8886E-02  5.5950E-02 -1.3844E-01  8.9477E-05
             5.0314E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2358
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.8705E-05  1.8472E-02 -3.1148E-02 -1.9059E-02 -3.3209E-02
 SE:             2.9859E-02  1.2930E-02  1.4672E-02  2.6813E-02  2.2063E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9843E-01  1.5312E-01  3.3754E-02  4.7722E-01  1.3227E-01

 ETASHRINKSD(%)  1.0000E-10  5.6682E+01  5.0847E+01  1.0172E+01  2.6087E+01
 ETASHRINKVR(%)  1.0000E-10  8.1235E+01  7.5840E+01  1.9308E+01  4.5369E+01
 EBVSHRINKSD(%)  4.1594E-01  6.8737E+01  5.4861E+01  6.9387E+00  1.9956E+01
 EBVSHRINKVR(%)  8.3015E-01  9.0226E+01  7.9624E+01  1.3396E+01  3.5930E+01
 RELATIVEINF(%)  9.8855E+01  2.2001E+00  5.0103E+00  2.0232E+01  1.5549E+01
 EPSSHRINKSD(%)  4.4577E+01
 EPSSHRINKVR(%)  6.9283E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1714.8819708719141     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -979.73114430817589     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    33.28
 Elapsed covariance  time in seconds:     7.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1714.882       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.27E-01  1.75E+00  1.57E+00  9.49E-01  9.58E-01  3.79E+00  1.22E+00  7.03E-01  1.22E+00  9.71E-01
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.02E-02  8.85E-02  9.88E-01  1.01E-01  2.71E-01  5.97E-02  8.76E-01  8.27E-01  6.92E-02  2.87E-01  6.20E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.14E-04
 
 TH 2
+       -2.61E-04  7.84E-03
 
 TH 3
+        2.97E-03 -4.60E-02  9.77E-01
 
 TH 4
+        4.78E-04 -7.17E-03  8.30E-02  1.01E-02
 
 TH 5
+        6.94E-04 -1.19E-02  2.65E-01  2.24E-02  7.36E-02
 
 TH 6
+        1.45E-04 -5.84E-05 -2.38E-03 -3.02E-04 -2.39E-04  3.56E-03
 
 TH 7
+        3.36E-03 -6.93E-02  4.06E-01  6.81E-02  1.08E-01  1.40E-03  7.67E-01
 
 TH 8
+        2.34E-03 -3.92E-02  7.92E-01  6.91E-02  2.13E-01 -3.33E-03  3.46E-01  6.85E-01
 
 TH 9
+       -4.12E-04  4.25E-03 -4.47E-02 -5.19E-03 -1.16E-02  2.19E-04 -2.90E-02 -3.68E-02  4.79E-03
 
 TH10
+        5.86E-04 -8.49E-03  2.55E-01  1.99E-02  7.21E-02  6.51E-04  8.40E-02  2.01E-01 -1.00E-02  8.22E-02
 
 TH11
+        2.87E-04  7.98E-04 -2.36E-03 -7.66E-04 -1.27E-03 -2.41E-04 -9.10E-03 -1.42E-03 -2.14E-05 -2.81E-03  3.84E-03
 
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
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.02E-02
 
 TH 2
+       -9.73E-02  8.85E-02
 
 TH 3
+        9.95E-02 -5.25E-01  9.88E-01
 
 TH 4
+        1.57E-01 -8.03E-01  8.34E-01  1.01E-01
 
 TH 5
+        8.46E-02 -4.94E-01  9.86E-01  8.18E-01  2.71E-01
 
 TH 6
+        8.04E-02 -1.11E-02 -4.04E-02 -5.03E-02 -1.47E-02  5.97E-02
 
 TH 7
+        1.27E-01 -8.93E-01  4.69E-01  7.72E-01  4.56E-01  2.67E-02  8.76E-01
 
 TH 8
+        9.34E-02 -5.35E-01  9.69E-01  8.29E-01  9.50E-01 -6.74E-02  4.78E-01  8.27E-01
 
 TH 9
+       -1.97E-01  6.94E-01 -6.53E-01 -7.45E-01 -6.15E-01  5.31E-02 -4.79E-01 -6.43E-01  6.92E-02
 
 TH10
+        6.76E-02 -3.35E-01  9.00E-01  6.87E-01  9.27E-01  3.80E-02  3.34E-01  8.48E-01 -5.07E-01  2.87E-01
 
 TH11
+        1.53E-01  1.45E-01 -3.86E-02 -1.23E-01 -7.55E-02 -6.51E-02 -1.68E-01 -2.77E-02 -5.00E-03 -1.58E-01  6.20E-02
 
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
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.25E+03
 
 TH 2
+       -2.67E+02  1.44E+03
 
 TH 3
+       -1.23E+01  3.19E+01  6.79E+01
 
 TH 4
+       -9.79E+01  3.71E+01 -1.19E+01  1.13E+03
 
 TH 5
+        6.16E+01 -1.72E+01 -1.68E+02 -1.54E+02  7.68E+02
 
 TH 6
+       -7.39E+01  7.41E+01  9.35E+00  6.23E+01 -3.88E+01  3.05E+02
 
 TH 7
+       -1.87E+01  1.01E+02  2.88E+00 -5.16E+01  1.81E+00  2.85E-01  1.17E+01
 
 TH 8
+        5.16E+00 -1.81E+00 -2.25E+01 -1.76E+01  4.63E+00  6.16E+00  2.22E-01  2.53E+01
 
 TH 9
+        2.01E+02 -6.78E+02  2.36E+01  3.41E+02 -1.16E+02 -2.28E+01 -5.86E+01 -5.62E+00  8.58E+02
 
 TH10
+       -3.38E-01 -1.28E+02 -2.63E+00  4.31E+01 -1.44E+02 -2.23E+01 -7.35E+00  6.85E+00  5.61E+01  1.22E+02
 
 TH11
+       -9.06E+01 -1.12E+02 -2.46E+01  7.49E+01  1.67E+01  1.07E+00 -5.31E+00 -5.76E-01  7.32E+01  5.97E+01  3.27E+02
 
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
 #CPUT: Total CPU Time in Seconds,       40.647
Stop Time:
Wed Sep 29 18:07:25 CDT 2021
