Thu Sep 30 10:05:11 CDT 2021
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
$DATA ../../../../data/spa2/D/dat93.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m93.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   37783.2976457204        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.6336E+02  8.2829E+02  4.8762E+01  7.5067E+02  1.1469E+02 -3.3243E+03 -1.5301E+03 -8.7076E+01 -2.2124E+03 -9.7726E+02
            -7.1427E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -349.122215983537        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2628E+00  1.2118E+00  9.6087E-01  1.3878E+00  1.0819E+00  1.8704E+00  1.1746E+00  9.7966E-01  1.0179E+00  9.6289E-01
             1.4802E+01
 PARAMETER:  3.3335E-01  2.9207E-01  6.0079E-02  4.2770E-01  1.7876E-01  7.2615E-01  2.6089E-01  7.9448E-02  1.1772E-01  6.2185E-02
             2.7948E+00
 GRADIENT:   1.6437E+01  2.9702E+01  1.0491E+00  5.8971E+01 -4.1886E+00  2.6615E+01 -1.8169E+01  3.2145E+00 -1.1300E+01  8.7888E+00
            -7.4665E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -482.135041849315        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.2704E+00  9.1264E-01  8.7434E+00  1.8936E+00  2.7592E+00  1.7780E+00  6.7694E+00  3.4569E-01  1.1635E+00  6.6268E-02
             1.5646E+01
 PARAMETER:  3.3937E-01  8.5843E-03  2.2683E+00  7.3849E-01  1.1149E+00  6.7547E-01  2.0124E+00 -9.6222E-01  2.5147E-01 -2.6140E+00
             2.8502E+00
 GRADIENT:   4.0304E+01  7.2785E+00 -1.4294E+00  1.9240E+01  1.8650E+00  1.1943E+01  4.3324E+01  1.2181E-02  1.2283E+01  1.2257E-02
             1.0474E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -496.443214357009        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      315
 NPARAMETR:  1.0764E+00  6.0120E-01  1.9392E+01  1.8489E+00  2.5675E+00  1.6042E+00  8.3924E+00  2.1391E-01  8.6142E-01  4.1112E-01
             1.4291E+01
 PARAMETER:  1.7364E-01 -4.0882E-01  3.0649E+00  7.1460E-01  1.0429E+00  5.7261E-01  2.2273E+00 -1.4422E+00 -4.9168E-02 -7.8886E-01
             2.7596E+00
 GRADIENT:  -8.1029E+00  7.4696E+00  2.6609E-01  2.0920E+01 -1.0448E+01  1.5581E+01  8.3080E+00  9.8610E-04 -4.6215E+00  5.6184E-01
             8.0974E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -498.554876604243        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      491
 NPARAMETR:  1.0765E+00  4.7738E-01  3.6408E+01  1.8622E+00  3.0279E+00  1.4880E+00  8.0744E+00  2.7261E-01  9.0652E-01  3.4963E-01
             1.4226E+01
 PARAMETER:  1.7373E-01 -6.3944E-01  3.6948E+00  7.2177E-01  1.2079E+00  4.9744E-01  2.1887E+00 -1.1997E+00  1.8550E-03 -9.5087E-01
             2.7551E+00
 GRADIENT:  -6.6585E+00  1.9784E+00 -1.8316E-02  8.2468E+00 -2.9159E-01 -2.4630E+00 -5.2523E+00  3.2521E-04 -2.9550E-01  3.2739E-01
             1.9531E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -501.010061761998        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      653
 NPARAMETR:  1.0811E+00  2.7531E-01  3.1685E+02  1.9491E+00  3.0753E+00  1.4926E+00  9.9400E+00  1.4652E+00  9.8570E-01  6.2639E-02
             1.4132E+01
 PARAMETER:  1.7800E-01 -1.1899E+00  5.8584E+00  7.6735E-01  1.2234E+00  5.0049E-01  2.3966E+00  4.8199E-01  8.5597E-02 -2.6704E+00
             2.7485E+00
 GRADIENT:   6.7856E+00  4.7556E+00  3.8378E-03  9.5118E+00  1.0189E+00  5.7994E+00  1.5629E+02  1.2993E-04  1.7165E+00  1.1255E-02
             3.0900E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -501.223395424527        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:      845             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0810E+00  2.4266E-01  2.5802E+02  1.9639E+00  3.0507E+00  1.4754E+00  1.0220E+01  1.4905E+00  9.7541E-01  2.4961E-02
             1.4188E+01
 PARAMETER:  1.7789E-01 -1.3161E+00  5.6530E+00  7.7492E-01  1.2154E+00  4.8893E-01  2.4244E+00  4.9909E-01  7.5103E-02 -3.5904E+00
             2.7524E+00
 GRADIENT:   6.2791E+00  4.5326E+00  6.3179E-03  1.1193E+01  5.6159E-01  4.5338E+00  1.6591E+02  2.1122E-04  7.4636E-01  1.8280E-03
             3.2468E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -501.397721680542        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1027
 NPARAMETR:  1.0849E+00  1.9102E-01  1.4893E+02  1.9917E+00  3.0584E+00  1.4662E+00  1.0687E+01  1.4560E+00  9.7526E-01  1.0000E-02
             1.4189E+01
 PARAMETER:  1.8146E-01 -1.5554E+00  5.1035E+00  7.8899E-01  1.2179E+00  4.8271E-01  2.4691E+00  4.7566E-01  7.4946E-02 -4.9635E+00
             2.7524E+00
 GRADIENT:   3.0081E+00  2.2455E-01  8.7776E-03 -3.1671E+00  3.6015E-02  4.3006E-02  2.3073E+01  6.3813E-04 -6.0797E-01  0.0000E+00
            -5.1910E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -501.415706774700        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1212
 NPARAMETR:  1.0840E+00  1.9716E-01  6.8784E+01  1.9973E+00  3.0609E+00  1.4685E+00  1.0603E+01  1.4533E+00  9.7908E-01  1.0000E-02
             1.4250E+01
 PARAMETER:  1.8066E-01 -1.5237E+00  4.3310E+00  7.9181E-01  1.2187E+00  4.8422E-01  2.4611E+00  4.7380E-01  7.8856E-02 -4.8856E+00
             2.7567E+00
 GRADIENT:   6.3693E-01  2.8166E-01  3.0132E-03 -6.1912E-01  1.5464E-01  1.1834E-02  2.1539E+01  3.0339E-03 -2.2225E-01  0.0000E+00
            -1.0924E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -501.425855306568        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1398
 NPARAMETR:  1.0835E+00  1.9170E-01  5.7492E+01  1.9998E+00  3.0387E+00  1.4689E+00  1.0685E+01  1.3138E+00  9.7972E-01  1.0000E-02
             1.4254E+01
 PARAMETER:  1.8017E-01 -1.5518E+00  4.1516E+00  7.9303E-01  1.2114E+00  4.8454E-01  2.4689E+00  3.7291E-01  7.9511E-02 -4.8856E+00
             2.7570E+00
 GRADIENT:   3.0723E-01  2.6214E-01  2.0084E-03 -3.2405E-01 -1.0623E-01  2.0610E-01  2.2337E+01  3.6535E-03 -3.1197E-01  0.0000E+00
            -9.6604E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -501.429521878807        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1582
 NPARAMETR:  1.0834E+00  1.8804E-01  5.7204E+01  2.0020E+00  3.0359E+00  1.4686E+00  1.0737E+01  1.1876E+00  9.7879E-01  1.0000E-02
             1.4254E+01
 PARAMETER:  1.8011E-01 -1.5711E+00  4.1466E+00  7.9417E-01  1.2105E+00  4.8428E-01  2.4737E+00  2.7195E-01  7.8559E-02 -4.8856E+00
             2.7570E+00
 GRADIENT:   2.7458E-01  2.3051E-01  3.0966E-03  2.0397E-01 -1.8663E-01  2.5538E-01  2.2750E+01  3.0361E-03 -4.9255E-01  0.0000E+00
            -1.3147E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -501.433402600483        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1769
 NPARAMETR:  1.0833E+00  1.8539E-01  5.4746E+01  2.0029E+00  3.0393E+00  1.4680E+00  1.0762E+01  1.0872E+00  9.8260E-01  1.0000E-02
             1.4259E+01
 PARAMETER:  1.8000E-01 -1.5853E+00  4.1027E+00  7.9461E-01  1.2116E+00  4.8387E-01  2.4760E+00  1.8360E-01  8.2446E-02 -4.8856E+00
             2.7574E+00
 GRADIENT:   1.3398E-01  1.5059E-01  8.9673E-04 -3.7423E-01 -5.9816E-02  1.2779E-01  2.2845E+01  2.7860E-03 -1.3871E-01  0.0000E+00
            -5.4608E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -501.435012541282        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1956             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0831E+00  1.8276E-01  5.3756E+01  2.0039E+00  3.0417E+00  1.4674E+00  1.0795E+01  1.0687E+00  9.8448E-01  1.0000E-02
             1.4262E+01
 PARAMETER:  1.7987E-01 -1.5996E+00  4.0845E+00  7.9508E-01  1.2124E+00  4.8349E-01  2.4790E+00  1.6649E-01  8.4357E-02 -4.8856E+00
             2.7576E+00
 GRADIENT:   6.0289E+00  3.9847E+00  1.0817E-03  1.8348E+01  3.7836E-01  4.2576E+00  1.8689E+02  2.8052E-03  1.9906E-01  0.0000E+00
             3.2894E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -501.435646461996        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2139
 NPARAMETR:  1.0832E+00  1.8098E-01  5.4109E+01  2.0048E+00  3.0435E+00  1.4672E+00  1.0823E+01  1.0308E+00  9.8574E-01  1.0000E-02
             1.4264E+01
 PARAMETER:  1.7989E-01 -1.6094E+00  4.0910E+00  7.9555E-01  1.2130E+00  4.8332E-01  2.4817E+00  1.3035E-01  8.5636E-02 -4.8856E+00
             2.7577E+00
 GRADIENT:   3.0780E-02  8.1377E-02 -2.8974E-04 -8.8201E-01  4.1768E-02  6.4171E-02  2.3348E+01  2.5632E-03  1.1776E-01  0.0000E+00
             4.0513E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -501.436078386898        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2327             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0832E+00  1.7972E-01  5.4591E+01  2.0063E+00  3.0414E+00  1.4669E+00  1.0842E+01  1.0188E+00  9.8434E-01  1.0000E-02
             1.4263E+01
 PARAMETER:  1.7993E-01 -1.6164E+00  4.0999E+00  7.9627E-01  1.2123E+00  4.8318E-01  2.4835E+00  1.1862E-01  8.4211E-02 -4.8856E+00
             2.7576E+00
 GRADIENT:   6.0314E+00  3.9967E+00  1.6518E-03  1.8926E+01  3.2805E-01  4.2410E+00  1.8858E+02  2.4757E-03  8.8425E-02  0.0000E+00
             3.2605E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -501.436143855119        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2508
 NPARAMETR:  1.0832E+00  1.7939E-01  5.4498E+01  2.0063E+00  3.0424E+00  1.4669E+00  1.0848E+01  1.0170E+00  9.8495E-01  1.0000E-02
             1.4264E+01
 PARAMETER:  1.7994E-01 -1.6182E+00  4.0982E+00  7.9631E-01  1.2126E+00  4.8313E-01  2.4839E+00  1.1684E-01  8.4835E-02 -4.8856E+00
             2.7577E+00
 GRADIENT:   3.5285E-02  7.0328E-02 -7.2354E-06 -3.1055E-01 -1.7008E-02  4.7755E-02  2.3508E+01  2.4706E-03 -2.2344E-02  0.0000E+00
            -3.2702E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -501.436170800502        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:     2655
 NPARAMETR:  1.0832E+00  1.7908E-01  5.4578E+01  2.0064E+00  3.0431E+00  1.4668E+00  1.0851E+01  1.0072E+00  9.8514E-01  1.0000E-02
             1.4265E+01
 PARAMETER:  1.7995E-01 -1.6199E+00  4.0996E+00  7.9635E-01  1.2129E+00  4.8309E-01  2.4843E+00  1.0715E-01  8.5027E-02 -4.8856E+00
             2.7578E+00
 GRADIENT:   6.0222E+00  3.9918E+00  1.4140E-03  1.8683E+01  3.6767E-01  4.2343E+00  1.8894E+02  2.4159E-03  1.7474E-01  0.0000E+00
             3.2858E+01

0ITERATION NO.:   82    OBJECTIVE VALUE:  -501.436170800502        NO. OF FUNC. EVALS.:  69
 CUMULATIVE NO. OF FUNC. EVALS.:     2724
 NPARAMETR:  1.0833E+00  1.7890E-01  5.4651E+01  2.0067E+00  3.0430E+00  1.4668E+00  1.0861E+01  1.0061E+00  9.8510E-01  1.0000E-02
             1.4265E+01
 PARAMETER:  1.7995E-01 -1.6199E+00  4.0996E+00  7.9635E-01  1.2129E+00  4.8309E-01  2.4843E+00  1.0715E-01  8.5027E-02 -4.8856E+00
             2.7578E+00
 GRADIENT:  -4.7858E-03  3.3380E-03 -1.9272E-05 -4.1184E-02  5.8627E-04 -2.6794E-04 -4.9961E-02  2.3409E-03  1.1451E-03  0.0000E+00
            -8.4619E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2724
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.7569E-02  7.8898E-02  2.2304E-04 -8.0576E-02  3.6683E-06
 SE:             2.7172E-02  1.9536E-02  6.0451E-05  1.4422E-02  2.6270E-05
 N:                     100         100         100         100         100

 P VAL.:         3.1029E-01  5.3817E-05  2.2473E-04  2.3139E-08  8.8895E-01

 ETASHRINKSD(%)  8.9718E+00  3.4551E+01  9.9797E+01  5.1686E+01  9.9912E+01
 ETASHRINKVR(%)  1.7139E+01  5.7164E+01  1.0000E+02  7.6657E+01  1.0000E+02
 EBVSHRINKSD(%)  1.5479E+01  3.6112E+01  9.9585E+01  4.4643E+01  9.9853E+01
 EBVSHRINKVR(%)  2.8562E+01  5.9183E+01  9.9998E+01  6.9356E+01  1.0000E+02
 RELATIVEINF(%)  6.5729E+01  2.8929E+01  2.9406E-04  1.3492E+01  3.3545E-05
 EPSSHRINKSD(%)  4.0340E+00
 EPSSHRINKVR(%)  7.9052E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -501.43617080050160     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       601.29006904510550     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    81.32
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    12.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -501.436       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.08E+00  1.79E-01  5.46E+01  2.01E+00  3.04E+00  1.47E+00  1.09E+01  1.01E+00  9.85E-01  1.00E-02  1.43E+01
 


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
+        4.87E+02
 
 TH 2
+       -6.34E+01  1.64E+02
 
 TH 3
+        3.69E-03 -9.08E-04  2.84E-07
 
 TH 4
+       -2.68E+02  3.38E+01 -5.80E-03  2.57E+02
 
 TH 5
+        1.01E+01 -2.60E+00  2.28E-04 -9.64E+00  3.74E-01
 
 TH 6
+       -1.26E+01 -2.27E+00  2.62E-03 -9.77E+00  5.07E-01  3.93E+01
 
 TH 7
+        4.85E+00  7.78E+00  1.40E-04 -5.82E+00  1.49E-01  4.67E-01  5.93E-01
 
 TH 8
+       -1.40E-01  6.60E-02 -2.97E-06  9.86E-02 -4.18E-03 -1.83E-02  4.08E-04  6.75E-05
 
 TH 9
+        3.91E+01 -9.21E+00  1.27E-03 -4.74E+01  1.82E+00  4.16E+00  9.11E-01 -1.82E-02  9.44E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.46E-01 -1.72E+00  3.64E-04 -7.40E+00  2.97E-01  2.92E+00  1.33E-01 -2.81E-03  1.87E+00  0.00E+00  6.24E-01
 
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
+        2.90E+02
 
 TH 2
+        3.66E+00  1.64E+02
 
 TH 3
+       -1.31E-03 -1.67E-04  7.25E-06
 
 TH 4
+       -3.91E+01  4.02E+00 -4.62E-03  1.40E+02
 
 TH 5
+        1.30E+00 -1.58E+00 -1.27E-03 -5.05E+00  2.38E+00
 
 TH 6
+       -1.08E+01  2.23E+00  4.51E-03 -1.66E+01  3.47E-01  5.29E+01
 
 TH 7
+        1.80E+00  8.78E+00  1.78E-04 -3.85E+00  1.47E-02  1.23E+00  8.35E-01
 
 TH 8
+       -7.16E-02  5.15E-02  1.53E-04  1.64E-02 -6.54E-03  7.48E-02  7.05E-04  5.22E-01
 
 TH 9
+       -3.34E+00 -6.19E+00 -2.03E-03 -2.86E+01  2.14E+00 -9.22E+00 -1.02E-01 -2.42E-01  4.31E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.93E+00 -2.06E+00  3.96E-04 -7.92E+00  2.36E-01  2.75E+00  2.05E-02 -2.10E-03  3.34E+00  0.00E+00  2.85E+00
 
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
+        3.24E+02
 
 TH 2
+        5.88E+01  1.75E+02
 
 TH 3
+       -9.96E-03  7.31E-04  2.69E-06
 
 TH 4
+        1.08E+02  8.05E+00 -2.29E-03  1.32E+02
 
 TH 5
+       -5.02E+00  2.91E+00 -7.43E-04 -7.24E+00  1.03E+00
 
 TH 6
+       -7.25E+00  3.78E+00 -1.98E-03 -2.14E+01  8.59E-01  4.46E+01
 
 TH 7
+       -1.27E+00  7.54E+00  8.62E-05 -3.38E+00  3.94E-01  4.44E-01  6.05E-01
 
 TH 8
+       -7.82E-04  1.36E-03  1.39E-07 -7.86E-04  5.69E-05  1.73E-04  7.22E-05  2.45E-07
 
 TH 9
+       -3.28E+01 -6.08E-01 -3.98E-04 -2.53E+01  2.00E+00  1.41E+01  1.06E-01  1.28E-04  2.65E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.54E+01 -2.20E+00 -9.38E-05 -2.15E+01  1.46E+00  6.87E+00  6.84E-01 -1.93E-06  5.27E+00  0.00E+00  3.14E+01
 
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
 #CPUT: Total CPU Time in Seconds,       87.622
Stop Time:
Thu Sep 30 10:06:47 CDT 2021
