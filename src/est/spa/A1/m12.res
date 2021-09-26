Sat Sep 25 07:53:04 CDT 2021
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
$DATA ../../../../data/spa/A1/dat12.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m12.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1388.16730593263        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.4390E+01 -4.3915E+01  5.0384E+00 -7.1736E+01  3.3944E+01  1.2561E+01 -3.2099E+01  8.7489E+00 -2.9736E+01 -3.2114E+01
            -4.5934E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1503.76275500073        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0036E+00  1.0567E+00  1.1063E+00  1.0115E+00  1.0514E+00  9.4932E-01  1.2103E+00  8.6129E-01  1.0446E+00  1.0581E+00
             1.7344E+00
 PARAMETER:  1.0361E-01  1.5515E-01  2.0098E-01  1.1141E-01  1.5014E-01  4.7992E-02  2.9083E-01 -4.9322E-02  1.4363E-01  1.5650E-01
             6.5066E-01
 GRADIENT:   3.5377E+01 -6.4535E+00  8.4136E+00 -2.3321E+01 -6.8074E+00 -3.7264E+00 -2.5982E+00  4.3691E+00  4.9841E+00 -5.3458E+00
            -7.3795E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1506.90517586795        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0036E+00  7.3056E-01  9.1987E-01  1.2685E+00  8.4812E-01  9.1891E-01  1.8122E+00  3.7204E-01  8.4813E-01  1.0846E+00
             1.6662E+00
 PARAMETER:  1.0361E-01 -2.1394E-01  1.6474E-02  3.3785E-01 -6.4737E-02  1.5438E-02  6.9457E-01 -8.8875E-01 -6.4723E-02  1.8118E-01
             6.1057E-01
 GRADIENT:   3.6077E+01  2.8495E+01 -3.5597E+01  9.0528E+01  2.6377E+01 -1.8781E+01  3.9111E+00  2.5980E+00 -3.1245E+00  1.7072E+01
            -1.3067E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1513.99180421220        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.8287E-01  5.7222E-01  7.3035E-01  1.2778E+00  6.4634E-01  9.6652E-01  2.0986E+00  9.0016E-02  8.1444E-01  7.3825E-01
             1.7008E+00
 PARAMETER:  8.2722E-02 -4.5824E-01 -2.1423E-01  3.4515E-01 -3.3643E-01  6.5944E-02  8.4126E-01 -2.3078E+00 -1.0525E-01 -2.0347E-01
             6.3107E-01
 GRADIENT:  -1.5127E+01  2.0656E+01  2.5133E+01  2.7977E+00 -3.4716E+01  1.2198E+00  5.0491E+00  1.8608E-01 -3.2854E-01 -2.8376E+00
            -2.5961E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1517.31504109472        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.8335E-01  3.1108E-01  6.8982E-01  1.4000E+00  5.6387E-01  9.6129E-01  2.8882E+00  1.0363E-02  7.9149E-01  7.5165E-01
             1.6914E+00
 PARAMETER:  8.3210E-02 -1.0677E+00 -2.7133E-01  4.3647E-01 -4.7293E-01  6.0519E-02  1.1606E+00 -4.4695E+00 -1.3384E-01 -1.8549E-01
             6.2558E-01
 GRADIENT:  -4.9379E+00  5.4802E+00  5.0655E+00  1.1349E+01 -1.2619E+01  2.5234E-01  4.2058E+00  2.9393E-03  5.0228E-01  1.3599E+00
            -1.0883E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1518.50830050420        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      438
 NPARAMETR:  9.8773E-01  2.3587E-01  8.3576E-01  1.4775E+00  6.3338E-01  9.6238E-01  3.3147E+00  1.0000E-02  7.8685E-01  8.5079E-01
             1.7031E+00
 PARAMETER:  8.7651E-02 -1.3445E+00 -7.9419E-02  4.9036E-01 -3.5669E-01  6.1656E-02  1.2984E+00 -5.0886E+00 -1.3972E-01 -6.1594E-02
             6.3247E-01
 GRADIENT:   3.0463E-01  2.8000E+00  3.6673E+00  7.4998E+00 -5.8680E+00  5.0601E-01  1.0842E+00  0.0000E+00 -3.8219E-01  2.6913E-01
            -6.8890E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1518.99480554437        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      615
 NPARAMETR:  9.8214E-01  1.0811E-01  8.5672E-01  1.5442E+00  6.1832E-01  9.6343E-01  4.3412E+00  1.0000E-02  7.9380E-01  8.8755E-01
             1.6976E+00
 PARAMETER:  8.1975E-02 -2.1246E+00 -5.4641E-02  5.3451E-01 -3.8075E-01  6.2740E-02  1.5681E+00 -8.2744E+00 -1.3093E-01 -1.9294E-02
             6.2921E-01
 GRADIENT:  -4.7778E+00 -2.0037E+00  7.0428E+00  5.1709E+00 -9.1047E+00  1.2333E+00 -6.9750E+00  0.0000E+00  7.2430E+00  1.3606E+00
            -1.0238E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1519.22604237629        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      790
 NPARAMETR:  9.8385E-01  8.3039E-02  8.7584E-01  1.5623E+00  6.2556E-01  9.5941E-01  4.8992E+00  1.0000E-02  7.7260E-01  9.0845E-01
             1.7058E+00
 PARAMETER:  8.3715E-02 -2.3884E+00 -3.2568E-02  5.4618E-01 -3.6911E-01  5.8566E-02  1.6891E+00 -9.3430E+00 -1.5800E-01  3.9876E-03
             6.3401E-01
 GRADIENT:   9.0165E-01  5.2121E-01 -3.3407E-01  1.0063E+00  7.1768E-01 -2.2921E-01  5.3096E-01  0.0000E+00 -7.8508E-01 -2.3076E-02
             4.1025E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1521.84152696047        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      970
 NPARAMETR:  9.7918E-01  2.7198E-02  8.4886E-01  1.5818E+00  5.9755E-01  9.6326E-01  7.7168E+00  1.0000E-02  7.4863E-01  9.1747E-01
             1.6863E+00
 PARAMETER:  7.8958E-02 -3.5046E+00 -6.3863E-02  5.5858E-01 -4.1491E-01  6.2564E-02  2.1434E+00 -1.4251E+01 -1.8951E-01  1.3870E-02
             6.2251E-01
 GRADIENT:  -7.7191E+00 -4.6426E+00  1.9000E+01  3.2035E+01 -3.3471E+01  1.6922E+00 -1.6048E+01  0.0000E+00 -2.5152E-01  1.5332E+01
            -1.5855E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1524.10975903557        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1148
 NPARAMETR:  9.7937E-01  1.0000E-02  8.0850E-01  1.5741E+00  5.7642E-01  9.6063E-01  1.1848E+01  1.0000E-02  7.3896E-01  8.9066E-01
             1.6959E+00
 PARAMETER:  7.9151E-02 -4.6140E+00 -1.1258E-01  5.5369E-01 -4.5092E-01  5.9838E-02  2.5722E+00 -1.9309E+01 -2.0251E-01 -1.5792E-02
             6.2822E-01
 GRADIENT:  -5.7919E+00  0.0000E+00  9.0019E+00  6.2844E+00 -1.9148E+01  7.9704E-01 -6.4142E+00  0.0000E+00 -6.8999E+00  1.2274E+01
             3.4463E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1524.58943808852        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1328
 NPARAMETR:  9.8202E-01  1.0000E-02  7.7790E-01  1.5603E+00  5.6561E-01  9.5777E-01  1.1675E+01  1.0000E-02  7.5595E-01  8.3591E-01
             1.7086E+00
 PARAMETER:  8.1861E-02 -4.6411E+00 -1.5116E-01  5.4490E-01 -4.6986E-01  5.6851E-02  2.5574E+00 -1.9567E+01 -1.7979E-01 -7.9240E-02
             6.3568E-01
 GRADIENT:   1.0364E+00  0.0000E+00 -2.7162E+00 -1.0325E+01  2.3225E+00 -4.0104E-01 -6.4701E+00  0.0000E+00 -5.2418E-01  5.0513E+00
             2.8418E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1524.63906211305        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1506
 NPARAMETR:  9.8184E-01  1.0000E-02  7.7888E-01  1.5622E+00  5.6559E-01  9.5803E-01  1.1683E+01  1.0000E-02  7.5884E-01  8.2933E-01
             1.7098E+00
 PARAMETER:  8.1675E-02 -4.6414E+00 -1.4990E-01  5.4607E-01 -4.6989E-01  5.7127E-02  2.5582E+00 -1.9557E+01 -1.7596E-01 -8.7138E-02
             6.3640E-01
 GRADIENT:   7.0855E+00  0.0000E+00 -2.1978E+01 -9.2739E+01  5.0320E+01  1.1950E+00  1.0760E+02  0.0000E+00 -5.1506E+01 -3.5805E+01
            -1.9904E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1524.67935302202        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1682
 NPARAMETR:  9.8150E-01  1.0000E-02  7.8239E-01  1.5681E+00  5.6642E-01  9.5876E-01  1.1722E+01  1.0000E-02  7.6428E-01  8.2302E-01
             1.7074E+00
 PARAMETER:  8.1323E-02 -4.6365E+00 -1.4540E-01  5.4986E-01 -4.6842E-01  5.7890E-02  2.5615E+00 -1.9499E+01 -1.6882E-01 -9.4780E-02
             6.3498E-01
 GRADIENT:  -2.0729E-01  0.0000E+00  5.7705E-01  2.6025E+00 -1.4553E+00 -5.4112E-02 -3.4178E+00  0.0000E+00  1.5540E+00  1.2212E+00
             7.0075E-01

0ITERATION NO.:   62    OBJECTIVE VALUE:  -1524.67935421536        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     1747
 NPARAMETR:  9.8158E-01  1.0000E-02  7.8229E-01  1.5684E+00  5.6631E-01  9.5886E-01  1.1711E+01  1.0000E-02  7.6434E-01  8.2295E-01
             1.7078E+00
 PARAMETER:  8.1327E-02 -4.6365E+00 -1.4539E-01  5.4985E-01 -4.6842E-01  5.7906E-02  2.5614E+00 -1.9499E+01 -1.6881E-01 -9.4880E-02
             6.3496E-01
 GRADIENT:  -1.3233E-01  0.0000E+00  3.5316E-01 -2.2403E+02  2.6810E+02 -3.0165E-02  4.2466E+01  0.0000E+00 -7.6056E+02 -7.7432E-01
            -2.0104E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1747
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4517E-04  1.4900E-02 -7.5656E-05 -1.0547E-02 -1.5433E-02
 SE:             2.9551E-02  7.9684E-03  1.9480E-04  2.7654E-02  2.2888E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9338E-01  6.1505E-02  6.9773E-01  7.0290E-01  5.0014E-01

 ETASHRINKSD(%)  1.0002E+00  7.3305E+01  9.9347E+01  7.3570E+00  2.3321E+01
 ETASHRINKVR(%)  1.9903E+00  9.2874E+01  9.9996E+01  1.4173E+01  4.1204E+01
 EBVSHRINKSD(%)  1.2296E+00  8.1599E+01  9.9281E+01  6.3863E+00  2.1216E+01
 EBVSHRINKVR(%)  2.4441E+00  9.6614E+01  9.9995E+01  1.2365E+01  3.7930E+01
 RELATIVEINF(%)  9.7413E+01  2.9773E+00  3.0604E-04  5.1241E+01  3.7216E+00
 EPSSHRINKSD(%)  3.7547E+01
 EPSSHRINKVR(%)  6.0996E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1524.6793542153605     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -789.52852765162231     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.33
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.83
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1524.679       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.82E-01  1.00E-02  7.82E-01  1.57E+00  5.66E-01  9.59E-01  1.17E+01  1.00E-02  7.64E-01  8.23E-01  1.71E+00
 


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
+        1.36E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -9.65E+02  0.00E+00  2.48E+06
 
 TH 4
+        1.48E+02  0.00E+00 -2.32E+03  4.18E+04
 
 TH 5
+       -5.66E+02  0.00E+00  5.62E+03  4.99E+03  4.45E+05
 
 TH 6
+        3.89E+01  0.00E+00 -1.76E+02  3.78E+01 -1.52E+02  2.23E+02
 
 TH 7
+       -4.66E+00  0.00E+00  6.20E+01  6.24E+01 -1.41E+02 -1.20E+00  3.44E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.30E+03  0.00E+00 -1.52E+04 -1.21E+03  4.13E+03  3.32E+02  3.12E+01  0.00E+00  1.93E+06
 
 TH10
+       -3.99E+06  0.00E+00  3.43E+06 -1.49E+03  4.93E+03 -4.08E+06  4.05E+01  0.00E+00 -1.09E+04  4.75E+06
 
 TH11
+        1.36E+02  0.00E+00 -1.81E+03 -3.89E+02  1.29E+03  4.11E+01  1.02E+01  0.00E+00 -1.05E+03 -1.23E+03  2.72E+04
 
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
 #CPUT: Total CPU Time in Seconds,       29.228
Stop Time:
Sat Sep 25 07:53:34 CDT 2021
