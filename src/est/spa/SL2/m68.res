Sat Sep 18 12:24:25 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat68.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m68.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1733.38426275356        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3250E+00  1.2642E+01 -2.6059E+01  5.9066E+01  3.0679E+01  3.7908E+00  2.2950E+01  9.2916E+00  4.6111E+01 -3.3844E+00
             6.1784E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1745.64288789069        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.8912E-01  1.0045E+00  1.0198E+00  1.0101E+00  9.9766E-01  9.8802E-01  8.6064E-01  9.4697E-01  8.0119E-01  1.0473E+00
             8.4929E-01
 PARAMETER:  8.9065E-02  1.0448E-01  1.1956E-01  1.1006E-01  9.7657E-02  8.7945E-02 -5.0077E-02  4.5510E-02 -1.2165E-01  1.4620E-01
            -6.3359E-02
 GRADIENT:  -1.9549E+00  6.1030E+01 -1.2494E+01  1.0123E+02  1.1718E+01 -4.8305E-01  4.9113E+00  4.9567E+00 -2.3290E-01 -2.6494E+00
             1.9458E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1747.03287724708        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      265
 NPARAMETR:  1.0107E+00  1.0753E+00  9.8493E-01  9.5996E-01  1.0248E+00  1.0129E+00  6.8924E-01  7.0893E-01  9.2067E-01  1.1516E+00
             8.4165E-01
 PARAMETER:  1.1060E-01  1.7257E-01  8.4812E-02  5.9134E-02  1.2451E-01  1.1286E-01 -2.7216E-01 -2.4399E-01  1.7349E-02  2.4119E-01
            -7.2393E-02
 GRADIENT:  -1.4747E+01  4.9791E+01 -4.7672E+00  7.7388E+01  5.3616E+00  4.0385E+00  1.7047E+00 -1.0613E+00  9.4904E+00  6.6301E+00
            -3.8382E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1749.41070080475        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  1.0181E+00  9.9898E-01  8.8822E-01  9.6622E-01  9.3144E-01  1.0028E+00  7.5460E-01  5.7761E-01  8.3026E-01  1.0257E+00
             8.4898E-01
 PARAMETER:  1.1794E-01  9.8975E-02 -1.8541E-02  6.5632E-02  2.8975E-02  1.0284E-01 -1.8157E-01 -4.4885E-01 -8.6019E-02  1.2536E-01
            -6.3714E-02
 GRADIENT:   1.7081E-01  3.0688E+00  1.8741E+00 -1.1031E+00 -2.0111E+00  6.8454E-02  6.2929E-02 -4.2959E-01 -2.5702E-01 -2.2791E-01
            -7.2541E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1749.87015831798        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      616
 NPARAMETR:  1.0148E+00  7.9816E-01  1.0561E+00  1.0991E+00  9.2895E-01  9.9943E-01  7.7402E-01  7.2229E-01  7.6683E-01  1.0522E+00
             8.5329E-01
 PARAMETER:  1.1468E-01 -1.2545E-01  1.5458E-01  1.9453E-01  2.6304E-02  9.9431E-02 -1.5616E-01 -2.2532E-01 -1.6549E-01  1.5090E-01
            -5.8661E-02
 GRADIENT:  -1.4107E+00  3.6678E+00  3.4043E-01  4.2616E+00  1.8169E+00 -8.4504E-02 -6.7598E-01 -6.8944E-01 -1.1692E+00 -1.2691E+00
             1.0559E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1750.68136907222        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      793
 NPARAMETR:  1.0100E+00  4.6528E-01  1.2696E+00  1.3144E+00  8.9732E-01  9.9448E-01  9.2045E-01  9.1894E-01  6.7549E-01  1.0663E+00
             8.5100E-01
 PARAMETER:  1.0994E-01 -6.6511E-01  3.3868E-01  3.7337E-01 -8.3435E-03  9.4468E-02  1.7102E-02  1.5466E-02 -2.9232E-01  1.6421E-01
            -6.1340E-02
 GRADIENT:  -5.9983E-01  5.8364E+00  2.0689E+00  1.8442E+01 -5.3023E+00  1.1432E-01 -7.1923E-02 -1.8928E-01  1.0882E+00 -5.1040E-01
            -7.5741E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1751.20513118886        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      970
 NPARAMETR:  1.0072E+00  2.9720E-01  1.4169E+00  1.4232E+00  9.0101E-01  9.9102E-01  1.1736E+00  1.0646E+00  6.2529E-01  1.0881E+00
             8.5077E-01
 PARAMETER:  1.0720E-01 -1.1133E+00  4.4847E-01  4.5292E-01 -4.2394E-03  9.0977E-02  2.6011E-01  1.6260E-01 -3.6954E-01  1.8444E-01
            -6.1615E-02
 GRADIENT:   1.3942E-01  4.1891E+00  1.5944E+00  2.1563E+01 -4.6730E+00  8.2300E-02 -5.7923E-02 -4.1454E-01 -8.4308E-01  7.3363E-01
            -1.7188E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1751.78651330357        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1146
 NPARAMETR:  1.0035E+00  1.2331E-01  1.6986E+00  1.5437E+00  9.3665E-01  9.8596E-01  2.0199E+00  1.3472E+00  5.8305E-01  1.1076E+00
             8.5203E-01
 PARAMETER:  1.0348E-01 -1.9930E+00  6.2981E-01  5.3419E-01  3.4550E-02  8.5856E-02  8.0307E-01  3.9802E-01 -4.3948E-01  2.0218E-01
            -6.0128E-02
 GRADIENT:  -1.4696E+00  2.2812E+00  1.5653E+00  3.0709E+01 -3.9592E+00 -4.9120E-01  1.8422E-02  1.1765E-01 -1.8375E-01 -7.7055E-01
             6.6691E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1752.19729865460        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1323
 NPARAMETR:  1.0021E+00  4.5736E-02  1.6062E+00  1.5808E+00  8.9436E-01  9.8604E-01  3.6309E+00  1.2749E+00  5.7052E-01  1.0959E+00
             8.5142E-01
 PARAMETER:  1.0213E-01 -2.9849E+00  5.7390E-01  5.5792E-01 -1.1646E-02  8.5937E-02  1.3895E+00  3.4286E-01 -4.6121E-01  1.9161E-01
            -6.0847E-02
 GRADIENT:  -1.0056E+00  4.3454E-01 -2.3024E+00  1.4559E+01 -1.3701E+00  4.9901E-02  1.7786E-02  6.9407E-01  5.5892E-01  1.7018E+00
             8.2248E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1752.37567720844        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1502
 NPARAMETR:  1.0021E+00  1.2854E-02  1.7812E+00  1.6075E+00  9.3203E-01  9.8507E-01  8.0292E+00  1.4177E+00  5.6192E-01  1.1113E+00
             8.4948E-01
 PARAMETER:  1.0214E-01 -4.2541E+00  6.7729E-01  5.7467E-01  2.9614E-02  8.4959E-02  2.1831E+00  4.4906E-01 -4.7639E-01  2.0549E-01
            -6.3132E-02
 GRADIENT:  -1.1314E-01  1.1630E-01  7.7276E-01  6.5199E+00 -8.8552E-01 -1.0641E-02  3.1242E-02 -9.2495E-02  8.9557E-02 -3.5874E-01
            -2.1353E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1752.39928074587        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1666
 NPARAMETR:  1.0020E+00  1.0000E-02  1.7657E+00  1.6068E+00  9.2863E-01  9.8493E-01  7.5845E+00  1.4066E+00  5.6156E-01  1.1112E+00
             8.4995E-01
 PARAMETER:  1.0199E-01 -4.5078E+00  6.6855E-01  5.7423E-01  2.5957E-02  8.4817E-02  2.1261E+00  4.4115E-01 -4.7703E-01  2.0540E-01
            -6.2576E-02
 GRADIENT:  -1.7891E-01  2.0860E-02  2.3309E-02 -6.0406E-01  4.6940E-01 -1.7313E-02  2.1394E-03  9.9157E-02 -1.7892E-01 -5.8731E-02
             1.0664E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1752.40301404568        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1842
 NPARAMETR:  1.0021E+00  1.0000E-02  1.7310E+00  1.6051E+00  9.1928E-01  9.8489E-01  4.1157E+00  1.3760E+00  5.6286E-01  1.1052E+00
             8.4964E-01
 PARAMETER:  1.0205E-01 -4.5315E+00  6.4871E-01  5.7321E-01  1.5833E-02  8.4773E-02  1.5148E+00  4.1921E-01 -4.7472E-01  2.0002E-01
            -6.2940E-02
 GRADIENT:   5.0814E-02  0.0000E+00 -2.4916E-02  3.5025E-01 -1.8581E-03 -5.2674E-02  1.3120E-03 -2.6973E-02  1.1336E-02 -8.0902E-03
            -7.9419E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1752.40426037456        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2020
 NPARAMETR:  1.0020E+00  1.0000E-02  1.7357E+00  1.6053E+00  9.2051E-01  9.8507E-01  8.4405E-01  1.3803E+00  5.6300E-01  1.1061E+00
             8.4980E-01
 PARAMETER:  1.0201E-01 -4.5359E+00  6.5143E-01  5.7329E-01  1.7173E-02  8.4957E-02 -6.9549E-02  4.2227E-01 -4.7447E-01  2.0084E-01
            -6.2752E-02
 GRADIENT:  -6.5922E-02  0.0000E+00  1.1158E-02 -1.2602E-02 -1.1125E-02  1.2930E-02  1.7033E-04  1.7290E-03 -8.1716E-03 -2.0221E-03
             7.1490E-03

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1752.40435185610        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2195
 NPARAMETR:  1.0020E+00  1.0000E-02  1.7351E+00  1.6052E+00  9.2036E-01  9.8511E-01  1.6082E-01  1.3797E+00  5.6303E-01  1.1060E+00
             8.4978E-01
 PARAMETER:  1.0203E-01 -4.5439E+00  6.5109E-01  5.7327E-01  1.7013E-02  8.5001E-02 -1.7275E+00  4.2190E-01 -4.7443E-01  2.0078E-01
            -6.2775E-02
 GRADIENT:   1.5528E-02  0.0000E+00  3.8576E-03 -3.2493E-02 -3.8905E-03  4.3496E-02  6.5780E-06  4.4479E-04 -5.5821E-03  1.7950E-03
            -2.5073E-03

0ITERATION NO.:   67    OBJECTIVE VALUE:  -1752.40435206568        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     2261
 NPARAMETR:  1.0020E+00  1.0000E-02  1.7351E+00  1.6051E+00  9.2036E-01  9.8500E-01  1.5986E-01  1.3797E+00  5.6305E-01  1.1060E+00
             8.4979E-01
 PARAMETER:  1.0202E-01 -4.5440E+00  6.5108E-01  5.7327E-01  1.7010E-02  8.4982E-02 -1.7317E+00  4.2189E-01 -4.7442E-01  2.0077E-01
            -6.2770E-02
 GRADIENT:  -6.7509E-03  0.0000E+00  3.2028E-03  2.7210E-01 -3.0643E-03  2.6029E-02  4.3740E-04  8.8178E-04 -5.5690E-03  2.2268E-03
            -3.5150E-04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2261
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.5755E-05 -6.4014E-05 -3.4597E-02 -8.8874E-03 -3.9941E-02
 SE:             2.9903E-02  3.5332E-05  1.7543E-02  2.9026E-02  2.1749E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9851E-01  7.0019E-02  4.8596E-02  7.5946E-01  6.6289E-02

 ETASHRINKSD(%)  1.0000E-10  9.9882E+01  4.1228E+01  2.7587E+00  2.7139E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  6.5458E+01  5.4413E+00  4.6912E+01
 EBVSHRINKSD(%)  2.9850E-01  9.9888E+01  4.5555E+01  3.1162E+00  2.2582E+01
 EBVSHRINKVR(%)  5.9610E-01  1.0000E+02  7.0358E+01  6.1352E+00  4.0065E+01
 RELATIVEINF(%)  9.7497E+01  5.2199E-06  7.0369E+00  4.2858E+00  1.0925E+01
 EPSSHRINKSD(%)  4.5626E+01
 EPSSHRINKVR(%)  7.0435E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1752.4043520656767     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1017.2535255019385     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.43
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.78
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1752.404       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.00E-02  1.74E+00  1.61E+00  9.20E-01  9.85E-01  1.60E-01  1.38E+00  5.63E-01  1.11E+00  8.50E-01
 


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
+        1.15E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.38E+00  0.00E+00  6.24E+01
 
 TH 4
+       -1.22E+01  0.00E+00 -3.96E+01  1.28E+03
 
 TH 5
+       -8.71E+00  0.00E+00 -1.53E+02 -8.47E+01  6.84E+02
 
 TH 6
+        2.37E+00  0.00E+00 -1.68E+00 -5.61E+00 -3.06E+01  1.95E+02
 
 TH 7
+       -6.04E+00  0.00E+00  5.58E-01 -1.17E+00  2.54E+00  1.17E+00  4.08E+00
 
 TH 8
+       -2.25E+00  0.00E+00 -2.01E+01 -4.95E+00 -6.84E+00  1.12E+00 -2.98E-01  2.44E+01
 
 TH 9
+        1.29E+01  0.00E+00  6.11E+00 -1.77E+00 -8.45E+00 -5.54E+00 -4.01E+00 -1.48E+00  5.64E+02
 
 TH10
+        5.08E+00  0.00E+00 -2.86E+00 -3.57E+00 -8.90E+01 -5.79E+00  2.89E+00  9.40E+00 -2.51E+00  6.91E+01
 
 TH11
+        1.59E+01  0.00E+00 -2.61E+00 -2.18E+01  3.32E+01  6.82E+00 -9.87E+00  6.15E+00  2.96E+01  1.36E+01  3.32E+02
 
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
 #CPUT: Total CPU Time in Seconds,       33.279
Stop Time:
Sat Sep 18 12:25:00 CDT 2021
