Thu Sep 30 10:12:07 CDT 2021
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
$DATA ../../../../data/spa2/D/dat99.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m99.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   38484.6699013301        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.1181E+02  7.9174E+02  4.5709E+01  7.8493E+02  8.8530E+01 -3.9608E+03 -1.7995E+03 -6.8870E+01 -2.3767E+03 -1.0131E+03
            -7.1759E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -317.065509566407        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2869E+00  1.1758E+00  9.5368E-01  1.2865E+00  1.1662E+00  1.7493E+00  1.2813E+00  9.8305E-01  1.0660E+00  9.7018E-01
             1.4797E+01
 PARAMETER:  3.5225E-01  2.6197E-01  5.2569E-02  3.5189E-01  2.5373E-01  6.5919E-01  3.4789E-01  8.2903E-02  1.6389E-01  6.9722E-02
             2.7945E+00
 GRADIENT:   3.5744E+01  2.8206E+00 -4.2929E+00  3.2058E+01  2.9568E+00  2.0352E+01 -2.3189E+01  3.1509E+00 -1.0320E+01  7.9999E+00
            -9.8013E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -460.252630801717        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.2121E+00  7.1867E-01  8.3053E+00  1.7725E+00  2.0358E+00  2.3612E+00  6.4435E+00  3.6283E-01  1.3290E+00  1.2916E-01
             1.4493E+01
 PARAMETER:  2.9235E-01 -2.3036E-01  2.2169E+00  6.7236E-01  8.1090E-01  9.5915E-01  1.9631E+00 -9.1383E-01  3.8445E-01 -1.9467E+00
             2.7737E+00
 GRADIENT:   1.1165E+01 -2.2475E+00  8.6740E-01  4.3158E+01 -1.6904E+01  6.4256E+01  3.2639E+01  1.4794E-02  6.4261E+00  8.0081E-02
             6.3994E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -472.916350033588        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      249
 NPARAMETR:  1.2012E+00  9.2898E-01  1.1598E+01  1.5288E+00  2.8037E+00  1.9424E+00  6.2089E+00  3.9345E-01  7.7131E-01  5.8582E-01
             1.4232E+01
 PARAMETER:  2.8333E-01  2.6331E-02  2.5508E+00  5.2446E-01  1.1309E+00  7.6391E-01  1.9260E+00 -8.3281E-01 -1.5967E-01 -4.3475E-01
             2.7555E+00
 GRADIENT:   9.5771E+00 -2.5681E+00 -4.2181E-01  1.9204E+01  2.2378E+00  1.7920E+01 -3.5145E+01  4.9792E-03  3.2721E-01  1.0474E+00
             9.3426E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -497.547831588844        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      427
 NPARAMETR:  1.0043E+00  7.4042E-01  2.7362E+01  1.5390E+00  3.3810E+00  1.1636E+00  8.0429E+00  5.3388E+00  4.7689E-01  4.3016E-01
             1.4763E+01
 PARAMETER:  1.0429E-01 -2.0054E-01  3.4091E+00  5.3114E-01  1.3182E+00  2.5154E-01  2.1848E+00  1.7750E+00 -6.4047E-01 -7.4360E-01
             2.7921E+00
 GRADIENT:  -1.0173E+01 -1.8650E+00 -3.2658E-01  4.0863E+00  1.5008E+00  1.0112E+01  9.9142E-01  9.7684E-02  5.1246E-02  3.8632E-01
             1.8071E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -497.697481871959        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  1.0128E+00  7.9461E-01  2.3057E+01  1.5012E+00  3.2072E+00  1.1255E+00  7.8439E+00  4.2733E+00  4.6411E-01  5.0595E-01
             1.4646E+01
 PARAMETER:  1.1270E-01 -1.2991E-01  3.2380E+00  5.0627E-01  1.2654E+00  2.1823E-01  2.1597E+00  1.5524E+00 -6.6764E-01 -5.8132E-01
             2.7842E+00
 GRADIENT:   3.3953E+00 -1.4881E+00 -3.3431E-01 -3.7910E+00  3.3802E-01  1.7189E+00  7.2734E-01  1.0053E-01  7.4370E-01  6.1628E-01
             7.9110E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -497.810721294905        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      785
 NPARAMETR:  1.0075E+00  8.9655E-01  2.2319E+01  1.4801E+00  3.1516E+00  1.1229E+00  7.6888E+00  9.5468E-01  4.3866E-01  3.3290E-01
             1.4603E+01
 PARAMETER:  1.0749E-01 -9.2023E-03  3.2055E+00  4.9213E-01  1.2479E+00  2.1595E-01  2.1398E+00  5.3618E-02 -7.2402E-01 -9.9990E-01
             2.7812E+00
 GRADIENT:  -1.0479E+00  7.7091E-01 -2.3336E-01  1.6748E-01 -2.2632E-02  4.1939E-01  6.9801E-01  5.7356E-03  5.7762E-01  2.7851E-01
             4.4029E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -497.958093756364        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      963
 NPARAMETR:  1.0053E+00  9.0798E-01  4.3211E+01  1.4528E+00  3.2100E+00  1.1164E+00  7.5987E+00  2.3927E-02  3.1606E-01  4.8072E-02
             1.4555E+01
 PARAMETER:  1.0532E-01  3.4623E-03  3.8661E+00  4.7349E-01  1.2663E+00  2.1014E-01  2.1280E+00 -3.6328E+00 -1.0518E+00 -2.9351E+00
             2.7779E+00
 GRADIENT:  -4.1849E-01  1.1758E-02 -9.8281E-02  8.5617E-01  3.0437E-01 -6.7070E-02 -1.9097E-01  7.7334E-07 -2.9407E-01  5.7625E-03
            -7.2975E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -498.067517820335        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1142
 NPARAMETR:  1.0053E+00  8.7927E-01  4.9111E+01  1.4630E+00  3.1936E+00  1.1161E+00  7.6518E+00  2.5917E-02  3.6457E-01  3.5488E-02
             1.4533E+01
 PARAMETER:  1.0533E-01 -2.8664E-02  3.9941E+00  4.8050E-01  1.2611E+00  2.0988E-01  2.1349E+00 -3.5529E+00 -9.0903E-01 -3.2386E+00
             2.7764E+00
 GRADIENT:   8.1360E-02 -2.8415E-01 -8.2793E-02 -4.0166E-01  1.7913E-01 -2.4727E-01 -6.2718E-02  7.3076E-07 -9.8973E-02  3.2141E-03
            -4.9198E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -498.369416503105        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1308
 NPARAMETR:  1.0060E+00  8.8314E-01  9.7657E+02  1.4712E+00  3.2181E+00  1.1165E+00  8.0849E+00  1.0000E-02  3.9870E-01  1.0000E-02
             1.4535E+01
 PARAMETER:  1.0602E-01 -2.4277E-02  6.9840E+00  4.8606E-01  1.2688E+00  2.1021E-01  2.1900E+00 -8.2954E+00 -8.1956E-01 -6.7757E+00
             2.7766E+00
 GRADIENT:   5.2943E+00  2.0187E+00 -2.7889E-03 -3.7520E+00  8.0867E-01  4.8879E-01  1.1461E+02  0.0000E+00  6.4918E-01  0.0000E+00
             3.6967E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -498.492490958220        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1495
 NPARAMETR:  1.0041E+00  8.1268E-01  1.8395E+03  1.4864E+00  3.2175E+00  1.1167E+00  8.0223E+00  1.0000E-02  3.8333E-01  1.0000E-02
             1.4523E+01
 PARAMETER:  1.0410E-01 -1.0742E-01  7.6173E+00  4.9638E-01  1.2686E+00  2.1036E-01  2.1822E+00 -8.2954E+00 -8.5886E-01 -6.7757E+00
             2.7757E+00
 GRADIENT:  -5.6666E-01 -1.2838E-01 -1.7257E-03 -1.1698E+00 -1.3305E-01  2.8396E-01  6.7854E+00  0.0000E+00 -2.9512E-01  0.0000E+00
            -6.9521E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -498.574915082001        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1660
 NPARAMETR:  1.0064E+00  8.1146E-01  3.7050E+03  1.4975E+00  3.2255E+00  1.1166E+00  8.1815E+00  1.0000E-02  4.1864E-01  1.0000E-02
             1.4533E+01
 PARAMETER:  1.0636E-01 -1.0892E-01  8.3174E+00  5.0382E-01  1.2711E+00  2.1026E-01  2.2019E+00 -8.2954E+00 -7.7074E-01 -6.7757E+00
             2.7764E+00
 GRADIENT:   1.2007E+00  6.8212E-01 -8.9003E-04 -3.3185E+00  6.6231E-02 -6.9588E-02  1.0627E+01  0.0000E+00  2.4425E-02  0.0000E+00
             9.9653E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -498.638929613913        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1849             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0059E+00  7.7169E-01  3.5323E+05  1.5155E+00  3.2180E+00  1.1157E+00  8.3914E+00  1.0000E-02  4.4302E-01  1.0000E-02
             1.4522E+01
 PARAMETER:  1.0588E-01 -1.5918E-01  1.2875E+01  5.1576E-01  1.2687E+00  2.0950E-01  2.2272E+00 -8.2954E+00 -7.1413E-01 -6.7757E+00
             2.7757E+00
 GRADIENT:  -4.5512E+01 -2.2030E+01  1.5293E-07 -3.4690E+01  7.4936E+00 -1.3959E+00  1.4948E+02  0.0000E+00  1.0353E+01  0.0000E+00
             9.9420E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -498.641171400576        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:     1975
 NPARAMETR:  1.0052E+00  7.7007E-01  1.4941E+05  1.5157E+00  3.2320E+00  1.1164E+00  8.3366E+00  1.0000E-02  4.5810E-01  1.0000E-02
             1.4539E+01
 PARAMETER:  1.0516E-01 -1.6128E-01  1.2014E+01  5.1589E-01  1.2731E+00  2.1008E-01  2.2207E+00 -8.2954E+00 -6.8067E-01 -6.7757E+00
             2.7768E+00
 GRADIENT:  -1.1275E+02 -1.0121E+00  3.1025E-06  3.1241E+01 -5.8027E-01  5.0647E+00  1.2145E+02  0.0000E+00 -3.0966E+00  0.0000E+00
             7.6616E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -498.652787998131        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2151
 NPARAMETR:  1.0068E+00  7.6087E-01  1.6753E+05  1.5201E+00  3.2428E+00  1.1168E+00  8.3839E+00  1.0000E-02  4.5781E-01  1.0000E-02
             1.4541E+01
 PARAMETER:  1.0674E-01 -1.7329E-01  1.2129E+01  5.1880E-01  1.2765E+00  2.1045E-01  2.2263E+00 -8.2954E+00 -6.8129E-01 -6.7757E+00
             2.7770E+00
 GRADIENT:  -3.2093E+02 -1.2491E+01  5.0308E-06  4.2658E+01  5.2072E+00  6.7940E+01  1.2637E+02  0.0000E+00 -1.5116E+00  0.0000E+00
             1.9615E+02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -498.657973656959        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2337
 NPARAMETR:  1.0058E+00  7.5249E-01  9.5601E+04  1.5222E+00  3.2332E+00  1.1160E+00  8.3504E+00  1.0000E-02  4.5592E-01  1.0000E-02
             1.4540E+01
 PARAMETER:  1.0581E-01 -1.8437E-01  1.1568E+01  5.2013E-01  1.2735E+00  2.0973E-01  2.2223E+00 -8.2954E+00 -6.8544E-01 -6.7757E+00
             2.7769E+00
 GRADIENT:  -5.2382E+00 -1.0493E+00 -3.7620E-05 -1.9790E+00  2.4760E-01 -3.6087E-01  1.2203E+01  0.0000E+00  1.6731E-01  0.0000E+00
             5.8291E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -498.666667669606        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     2508
 NPARAMETR:  1.0055E+00  7.4003E-01  9.9282E+04  1.5264E+00  3.2231E+00  1.1165E+00  8.4135E+00  1.0000E-02  4.5849E-01  1.0000E-02
             1.4532E+01
 PARAMETER:  1.0549E-01 -2.0106E-01  1.1606E+01  5.2294E-01  1.2704E+00  2.1016E-01  2.2298E+00 -8.2954E+00 -6.7982E-01 -6.7757E+00
             2.7764E+00
 GRADIENT:  -1.9664E+01 -2.2746E+00 -3.2522E-05 -1.2637E+00  3.9753E-01  5.7743E-01  1.3873E+01  0.0000E+00  4.8050E-01  0.0000E+00
             1.2799E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -498.668760868004        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2695
 NPARAMETR:  1.0069E+00  7.3573E-01  8.8351E+04  1.5277E+00  3.2279E+00  1.1187E+00  8.4214E+00  1.0000E-02  4.6441E-01  1.0000E-02
             1.4542E+01
 PARAMETER:  1.0683E-01 -2.0689E-01  1.1489E+01  5.2376E-01  1.2718E+00  2.1217E-01  2.2308E+00 -8.2954E+00 -6.6700E-01 -6.7757E+00
             2.7771E+00
 GRADIENT:  -8.4784E-01 -1.0468E+00 -3.9975E-05 -2.8908E+00  1.8776E-01  1.0343E+00  1.3093E+01  0.0000E+00  1.9854E-01  0.0000E+00
             3.8692E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -498.671981381845        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     2866
 NPARAMETR:  1.0058E+00  7.2876E-01  1.3150E+05  1.5305E+00  3.2338E+00  1.1176E+00  8.4571E+00  1.0000E-02  4.7108E-01  1.0000E-02
             1.4542E+01
 PARAMETER:  1.0578E-01 -2.1641E-01  1.1887E+01  5.2560E-01  1.2737E+00  2.1115E-01  2.2350E+00 -8.2954E+00 -6.5272E-01 -6.7757E+00
             2.7771E+00
 GRADIENT:  -3.2364E+02 -1.7139E+01 -1.3945E-07  2.5775E+01  1.0213E+01  9.0881E+01  1.8603E+01  0.0000E+00  3.8006E+00  0.0000E+00
             1.8233E+02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -498.672984922091        NO. OF FUNC. EVALS.: 210
 CUMULATIVE NO. OF FUNC. EVALS.:     3076             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0061E+00  7.2561E-01  1.3825E+05  1.5325E+00  3.2288E+00  1.1165E+00  8.4637E+00  1.0000E-02  4.7023E-01  1.0000E-02
             1.4536E+01
 PARAMETER:  1.0609E-01 -2.2074E-01  1.1937E+01  5.2689E-01  1.2721E+00  2.1022E-01  2.2358E+00 -8.2954E+00 -6.5453E-01 -6.7757E+00
             2.7766E+00
 GRADIENT:  -9.0736E+02 -3.1190E+01  4.2133E-05  1.3406E+02  2.4371E+01  2.7177E+02  1.3156E+02  0.0000E+00  5.0807E+00  0.0000E+00
             4.7703E+02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -498.673253634152        NO. OF FUNC. EVALS.: 203
 CUMULATIVE NO. OF FUNC. EVALS.:     3279
 NPARAMETR:  1.0061E+00  7.2439E-01  1.0386E+05  1.5338E+00  3.2261E+00  1.1165E+00  8.4691E+00  1.0000E-02  4.7275E-01  1.0000E-02
             1.4534E+01
 PARAMETER:  1.0612E-01 -2.2243E-01  1.1651E+01  5.2776E-01  1.2713E+00  2.1021E-01  2.2364E+00 -8.2954E+00 -6.4918E-01 -6.7757E+00
             2.7765E+00
 GRADIENT:  -1.2581E+02 -8.1229E-01  2.6656E-05  3.8418E+01 -2.4029E+00  5.5760E+00  1.1085E+01  0.0000E+00 -5.3836E+00  0.0000E+00
             4.1992E+01

0ITERATION NO.:  103    OBJECTIVE VALUE:  -498.673722019553        NO. OF FUNC. EVALS.: 108
 CUMULATIVE NO. OF FUNC. EVALS.:     3387
 NPARAMETR:  1.0061E+00  7.2391E-01  8.8633E+04  1.5353E+00  3.2254E+00  1.1166E+00  8.4744E+00  1.0000E-02  4.7377E-01  1.0000E-02
             1.4535E+01
 PARAMETER:  1.0623E-01 -2.2301E-01  1.1472E+01  5.2739E-01  1.2727E+00  2.1026E-01  2.2371E+00 -8.2954E+00 -6.4063E-01 -6.7757E+00
             2.7768E+00
 GRADIENT:   7.4872E-02  3.8332E-03 -1.6462E-04 -6.4633E-01  3.9465E-02  2.1047E-03  4.6241E-03  0.0000E+00  6.9856E-02  0.0000E+00
             1.9255E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3387
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.2854E-02  4.5331E-02  1.1714E-09 -4.7574E-02  2.6187E-05
 SE:             2.4341E-02  2.2795E-02  2.8526E-10  8.2494E-03  2.4328E-05
 N:                     100         100         100         100         100

 P VAL.:         1.7710E-01  4.6740E-02  4.0237E-05  8.0935E-09  2.8175E-01

 ETASHRINKSD(%)  1.8454E+01  2.3635E+01  1.0000E+02  7.2363E+01  9.9918E+01
 ETASHRINKVR(%)  3.3502E+01  4.1683E+01  1.0000E+02  9.2362E+01  1.0000E+02
 EBVSHRINKSD(%)  2.3409E+01  1.6910E+01  1.0000E+02  7.5867E+01  9.9863E+01
 EBVSHRINKVR(%)  4.1338E+01  3.0960E+01  1.0000E+02  9.4176E+01  1.0000E+02
 RELATIVEINF(%)  5.6071E+01  4.3968E+01  0.0000E+00  3.3578E+00  1.3816E-04
 EPSSHRINKSD(%)  2.2700E+00
 EPSSHRINKVR(%)  4.4885E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -498.67372201955277     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       604.05251782605433     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    91.57
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    12.86
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -498.674       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  7.24E-01  8.69E+04  1.53E+00  3.23E+00  1.12E+00  8.47E+00  1.00E-02  4.77E-01  1.00E-02  1.45E+01
 


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
+        5.22E+02
 
 TH 2
+       -4.50E+00  2.74E+00
 
 TH 3
+       -2.47E-06 -5.40E-08  1.40E-14
 
 TH 4
+       -2.08E+02  2.20E+01  5.05E-07  2.63E+02
 
 TH 5
+        5.94E+00 -7.76E-01 -9.86E-09 -8.48E+00  2.79E-01
 
 TH 6
+       -7.34E-01  7.96E+00 -3.68E-07  6.53E+00 -8.65E-01  1.22E+02
 
 TH 7
+        7.13E+00 -6.45E-01 -1.98E-08 -8.04E+00  2.58E-01 -1.96E-01  2.47E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.11E+01 -6.29E+00 -3.83E-08 -6.41E+01  2.14E+00 -1.14E+01  1.94E+00  0.00E+00  1.67E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.00E+00 -1.28E+00  5.50E-08 -9.84E+00  3.45E-01 -2.81E-01  2.73E-01  0.00E+00  2.74E+00  0.00E+00  8.28E-01
 
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
+        4.79E+02
 
 TH 2
+        6.78E+00  2.63E+01
 
 TH 3
+       -2.49E-06 -2.11E-07  1.92E-13
 
 TH 4
+       -6.41E+01  1.87E+01 -2.26E-07  2.08E+02
 
 TH 5
+        1.35E+00 -7.74E-01  5.22E-08 -6.82E+00  2.22E+00
 
 TH 6
+       -3.37E+01  4.45E+00 -1.03E-07 -1.70E+00 -3.56E-01  9.37E+01
 
 TH 7
+        2.73E+00  2.36E+00 -2.46E-08 -6.52E+00  1.41E-01 -7.21E-02  1.46E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.83E+00 -4.20E+00 -1.95E-07 -4.98E+01  1.73E+00 -4.66E+00  9.24E-01  0.00E+00  2.83E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.28E+01 -1.44E+00  5.49E-09 -1.11E+01  3.52E-01  1.40E+00  2.27E-01  0.00E+00  3.20E+00  0.00E+00  3.39E+00
 
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
+        5.36E+02
 
 TH 2
+        3.18E+01  2.17E+01
 
 TH 3
+       -4.33E-09 -9.09E-10  9.97E-18
 
 TH 4
+        7.24E+01  2.13E+01 -6.61E-09  1.95E+02
 
 TH 5
+       -5.62E+00 -1.08E-01  1.47E-10 -7.16E+00  6.75E-01
 
 TH 6
+       -6.78E+01 -4.86E+00  5.61E-09 -1.65E+01  6.81E-01  7.53E+01
 
 TH 7
+       -2.68E+00  2.47E+00  2.25E-10 -6.67E+00  5.30E-01  1.18E+00  1.54E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -5.82E+00 -9.19E-01  3.84E-10 -3.87E+01  1.16E+00 -3.44E+00  7.97E-01  0.00E+00  1.61E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.34E+01 -6.28E+00  3.64E-09 -2.58E+01  1.26E+00  7.15E+00  6.62E-01  0.00E+00  4.47E+00  0.00E+00  2.99E+01
 
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
 #CPUT: Total CPU Time in Seconds,      104.504
Stop Time:
Thu Sep 30 10:13:53 CDT 2021
