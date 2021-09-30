Wed Sep 29 18:29:13 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat78.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1672.24140651003        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5796E+02 -7.7691E+01 -6.2133E+01 -2.4811E+01  4.0207E+01  3.2928E+01 -1.5768E+01  1.8324E+01  1.1521E+01  6.6149E+00
             2.0203E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1691.85099494130        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.9308E-01  1.3463E+00  1.7512E+00  8.7179E-01  1.3871E+00  9.9019E-01  1.2327E+00  7.6168E-01  8.9401E-01  1.0700E+00
             9.4587E-01
 PARAMETER:  9.3055E-02  3.9733E-01  6.6029E-01 -3.7204E-02  4.2720E-01  9.0140E-02  3.0922E-01 -1.7223E-01 -1.2042E-02  1.6767E-01
             4.4352E-02
 GRADIENT:   4.7454E+01  9.9717E+00  7.7600E+00 -1.0160E+01  2.6627E+01 -3.0074E+00  1.9076E+01 -1.4328E+00 -4.0227E+00 -4.7698E+01
            -2.2221E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1697.94323244186        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.6924E-01  1.3516E+00  2.5395E+00  8.4373E-01  1.5017E+00  1.0115E+00  9.6074E-01  5.7801E-01  1.1146E+00  1.3874E+00
             9.8671E-01
 PARAMETER:  6.8760E-02  4.0130E-01  1.0320E+00 -6.9926E-02  5.0662E-01  1.1147E-01  5.9948E-02 -4.4817E-01  2.0846E-01  4.2746E-01
             8.6624E-02
 GRADIENT:  -8.5566E+00 -3.8698E+01  3.7247E+00 -4.0528E+01  2.3724E+00  6.3810E+00  4.2299E+00 -1.5578E-01 -1.7590E+00  1.1046E+00
             4.1395E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1699.57406240351        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  9.7423E-01  1.4031E+00  1.7106E+00  8.4283E-01  1.4060E+00  9.9689E-01  9.2276E-01  2.8899E-01  1.1439E+00  1.3009E+00
             9.6936E-01
 PARAMETER:  7.3892E-02  4.3870E-01  6.3685E-01 -7.0990E-02  4.4077E-01  9.6883E-02  1.9618E-02 -1.1414E+00  2.3449E-01  3.6306E-01
             6.8882E-02
 GRADIENT:   2.7239E+00  5.5040E+00 -2.5727E+00  8.9218E+00  5.9657E+00  6.4606E-01 -4.6998E-01  6.5223E-02 -3.8375E-01 -2.0157E+00
            -3.9082E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1699.65994993470        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      701
 NPARAMETR:  9.7368E-01  1.4589E+00  1.8240E+00  8.0324E-01  1.4344E+00  9.9484E-01  8.7092E-01  3.0738E-01  1.2175E+00  1.3350E+00
             9.7338E-01
 PARAMETER:  7.3329E-02  4.7769E-01  7.0105E-01 -1.1910E-01  4.6073E-01  9.4826E-02 -3.8208E-02 -1.0797E+00  2.9680E-01  3.8890E-01
             7.3016E-02
 GRADIENT:   9.1645E-01  3.8380E+00  8.1661E-01  4.0667E+00 -2.6695E+00 -2.1939E-01  4.2321E-01  1.2812E-02  5.7540E-01  5.4557E-01
             8.4444E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1699.66009238786        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.7369E-01  1.4705E+00  1.8189E+00  7.9533E-01  1.4374E+00  9.9477E-01  8.6331E-01  3.0477E-01  1.2283E+00  1.3369E+00
             9.7351E-01
 PARAMETER:  7.3340E-02  4.8558E-01  6.9821E-01 -1.2900E-01  4.6287E-01  9.4752E-02 -4.6987E-02 -1.0882E+00  3.0567E-01  3.9033E-01
             7.3155E-02
 GRADIENT:   8.2517E-01  3.6533E+00  8.3061E-01  3.9585E+00 -2.7846E+00 -2.6744E-01  2.6845E-01  1.2342E-02  5.4119E-01  6.4099E-01
             8.3928E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1699.69154585392        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1067
 NPARAMETR:  9.7273E-01  1.4797E+00  1.7481E+00  7.8655E-01  1.4329E+00  9.9549E-01  8.7203E-01  1.0648E-01  1.2203E+00  1.3296E+00
             9.7253E-01
 PARAMETER:  7.2353E-02  4.9181E-01  6.5855E-01 -1.4010E-01  4.5970E-01  9.5477E-02 -3.6929E-02 -2.1398E+00  2.9912E-01  3.8488E-01
             7.2141E-02
 GRADIENT:  -1.4060E+00  4.5522E-01  3.9900E-01  1.7239E+00 -1.0114E+00  1.2038E-02  6.5184E-02  2.2520E-03 -6.0171E-02  1.9648E-01
            -7.9301E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1699.69423595915        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1243
 NPARAMETR:  9.7274E-01  1.4948E+00  1.6879E+00  7.7649E-01  1.4290E+00  9.9548E-01  8.8014E-01  3.8800E-02  1.2190E+00  1.3235E+00
             9.7150E-01
 PARAMETER:  7.2363E-02  5.0196E-01  6.2348E-01 -1.5297E-01  4.5697E-01  9.5475E-02 -2.7671E-02 -3.1493E+00  2.9800E-01  3.8030E-01
             7.1089E-02
 GRADIENT:  -1.4898E+00  5.3227E-01 -2.0240E-02  2.1038E+00 -6.8634E-01 -8.0057E-03  2.1106E-01  6.1875E-04 -6.6725E-02 -8.6287E-02
            -3.1293E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1699.69433115138        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1423
 NPARAMETR:  9.7282E-01  1.5029E+00  1.6713E+00  7.7092E-01  1.4292E+00  9.9551E-01  8.8043E-01  2.4473E-02  1.2218E+00  1.3228E+00
             9.7125E-01
 PARAMETER:  7.2448E-02  5.0739E-01  6.1358E-01 -1.6017E-01  4.5711E-01  9.5497E-02 -2.7342E-02 -3.6102E+00  3.0036E-01  3.7974E-01
             7.0830E-02
 GRADIENT:  -1.3669E+00  3.7846E-01 -8.6664E-02  2.0316E+00 -7.4604E-01 -9.4874E-03  2.5155E-01  2.8098E-04 -8.3685E-02 -1.1557E-01
            -4.8563E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1699.69435519105        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1614
 NPARAMETR:  9.7287E-01  1.5088E+00  1.6614E+00  7.6690E-01  1.4297E+00  9.9552E-01  8.8004E-01  1.7724E-02  1.2245E+00  1.3225E+00
             9.7112E-01
 PARAMETER:  7.2492E-02  5.1135E-01  6.0765E-01 -1.6540E-01  4.5749E-01  9.5514E-02 -2.7788E-02 -3.9329E+00  3.0254E-01  3.7956E-01
             7.0696E-02
 GRADIENT:  -1.3176E+00  3.1811E-01 -1.1873E-01  2.0159E+00 -7.8208E-01 -1.0945E-02  2.7135E-01  1.5783E-04 -9.6249E-02 -1.3686E-01
            -6.1019E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1699.70092927104        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1799
 NPARAMETR:  9.7343E-01  1.5111E+00  1.6670E+00  7.6439E-01  1.4319E+00  9.9554E-01  8.7474E-01  1.0000E-02  1.2293E+00  1.3240E+00
             9.7134E-01
 PARAMETER:  7.3071E-02  5.1282E-01  6.1103E-01 -1.6868E-01  4.5902E-01  9.5526E-02 -3.3823E-02 -4.8569E+00  3.0643E-01  3.8067E-01
             7.0917E-02
 GRADIENT:  -5.1766E-02 -1.0701E+00 -9.1142E-02  1.0109E+00 -6.0263E-01 -3.5310E-03  5.6804E-02  0.0000E+00 -2.3142E-01 -1.0541E-01
            -5.6923E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1699.70157304936        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1976
 NPARAMETR:  9.7298E-01  1.5232E+00  1.6770E+00  7.5761E-01  1.4380E+00  9.9545E-01  8.6607E-01  1.0000E-02  1.2433E+00  1.3277E+00
             9.7185E-01
 PARAMETER:  7.2606E-02  5.2079E-01  6.1702E-01 -1.7759E-01  4.6329E-01  9.5442E-02 -4.3787E-02 -6.0540E+00  3.1777E-01  3.8348E-01
             7.1441E-02
 GRADIENT:  -1.2236E+00  4.1900E-01  6.1378E-02  2.2215E+00 -5.5718E-01 -6.9963E-02 -2.3165E-03  0.0000E+00 -1.1595E-01 -9.6028E-02
            -5.4926E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1699.70816046814        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2163             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7403E-01  1.5275E+00  1.6775E+00  7.5151E-01  1.4411E+00  9.9568E-01  8.6262E-01  1.0000E-02  1.2504E+00  1.3297E+00
             9.7208E-01
 PARAMETER:  7.3685E-02  5.2363E-01  6.1728E-01 -1.8567E-01  4.6541E-01  9.5666E-02 -4.7786E-02 -6.8549E+00  3.2346E-01  3.8497E-01
             7.1678E-02
 GRADIENT:   4.3460E+02  3.9533E+02  1.0316E+00  6.6220E+01  2.0759E+01  3.8011E+01  5.1061E+00  0.0000E+00  1.5081E+01  4.7381E+00
             8.7974E-01

0ITERATION NO.:   62    OBJECTIVE VALUE:  -1699.70816046814        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     2228
 NPARAMETR:  9.7429E-01  1.5289E+00  1.6748E+00  7.5248E-01  1.4410E+00  9.9586E-01  8.6175E-01  1.0000E-02  1.2518E+00  1.3300E+00
             9.7189E-01
 PARAMETER:  7.3685E-02  5.2363E-01  6.1728E-01 -1.8567E-01  4.6541E-01  9.5666E-02 -4.7786E-02 -6.8549E+00  3.2346E-01  3.8497E-01
             7.1678E-02
 GRADIENT:  -5.0153E-01 -1.1179E+00  8.5882E-02 -4.1657E-01  3.0526E-02 -6.1054E-02  3.8847E-02  0.0000E+00 -1.3046E-01 -3.1742E-02
             8.0092E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2228
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.5417E-04 -1.9038E-02 -7.8681E-05  7.2116E-03 -3.7886E-02
 SE:             2.9768E-02  1.9667E-02  4.7923E-05  2.2206E-02  2.3926E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9051E-01  3.3302E-01  1.0063E-01  7.4536E-01  1.1332E-01

 ETASHRINKSD(%)  2.7359E-01  3.4114E+01  9.9839E+01  2.5607E+01  1.9845E+01
 ETASHRINKVR(%)  5.4643E-01  5.6591E+01  1.0000E+02  4.4657E+01  3.5752E+01
 EBVSHRINKSD(%)  4.2573E-01  3.2333E+01  9.9839E+01  2.7181E+01  1.6628E+01
 EBVSHRINKVR(%)  8.4965E-01  5.4211E+01  1.0000E+02  4.6974E+01  3.0491E+01
 RELATIVEINF(%)  9.8829E+01  9.4073E-01  6.1503E-05  1.1135E+00  1.8153E+01
 EPSSHRINKSD(%)  4.0944E+01
 EPSSHRINKVR(%)  6.5124E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1699.7081604681441     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -964.55733390440594     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.12
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.46
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1699.708       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.74E-01  1.53E+00  1.68E+00  7.52E-01  1.44E+00  9.96E-01  8.63E-01  1.00E-02  1.25E+00  1.33E+00  9.72E-01
 


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
+        1.17E+03
 
 TH 2
+       -7.97E+00  3.01E+02
 
 TH 3
+        4.16E-01  1.03E+01  1.15E+01
 
 TH 4
+       -5.91E+00  4.07E+02 -1.80E+01  6.71E+02
 
 TH 5
+       -2.39E+00 -5.80E+01 -3.36E+01  3.50E+01  2.12E+02
 
 TH 6
+        1.58E+00 -2.07E+00  1.76E-01 -3.07E+00 -1.81E+00  1.97E+02
 
 TH 7
+        1.61E+00 -5.57E+00  6.51E+00 -8.29E+00 -9.62E+00 -5.59E-02  6.09E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.09E+00 -7.84E+00 -3.01E+00  2.49E+01  3.13E+00 -3.29E-01  3.43E+01  0.00E+00  4.26E+01
 
 TH10
+       -1.41E-01 -3.30E+00 -2.84E+00 -3.28E+00 -3.30E+01  5.76E-01 -3.97E-01  0.00E+00  1.29E+00  5.97E+01
 
 TH11
+       -1.04E+01 -1.50E+01 -8.95E+00 -4.82E+00  6.87E+00  1.43E+00  5.92E+00  0.00E+00  3.49E+00  1.72E+01  2.60E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       36.656
Stop Time:
Wed Sep 29 18:29:51 CDT 2021
