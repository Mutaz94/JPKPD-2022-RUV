Wed Sep 29 23:50:31 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat96.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m96.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1295.34141785536        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9151E+02  2.4773E+01  4.4783E+01  3.3538E+01  8.7463E+01  3.3663E+01 -2.1295E+00 -7.5471E+01 -4.6595E+01 -2.9659E+01
            -1.3754E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1719.20289525528        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1001E+00  1.0999E+00  1.3470E+00  1.0061E+00  1.1145E+00  1.1656E+00  8.3502E-01  7.3983E-01  6.9423E-01  6.5278E-01
             2.5046E+00
 PARAMETER:  1.9539E-01  1.9521E-01  3.9790E-01  1.0607E-01  2.0845E-01  2.5324E-01 -8.0299E-02 -2.0133E-01 -2.6495E-01 -3.2652E-01
             1.0181E+00
 GRADIENT:   3.2371E+02  4.8275E+01  1.1422E+01  3.9579E+01 -2.7652E+01  4.9008E+01 -1.2419E+01  1.1559E+00 -2.6118E+01  4.4452E+00
             5.0104E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1739.53743182018        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0003E+00  1.1954E+00  9.5106E-01  9.4460E-01  1.0237E+00  9.7198E-01  8.6751E-01  1.5544E-01  1.0647E+00  5.0257E-01
             2.5421E+00
 PARAMETER:  1.0030E-01  2.7844E-01  4.9825E-02  4.3006E-02  1.2347E-01  7.1582E-02 -4.2127E-02 -1.7615E+00  1.6267E-01 -5.8802E-01
             1.0330E+00
 GRADIENT:   1.2327E+02  5.6257E+01  1.6273E+00  5.0833E+01 -9.3389E+00  7.3891E+00  1.2272E+01  1.4682E-01  2.1888E+01  2.9582E+00
             8.5137E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1746.27870923538        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  9.5232E-01  1.0095E+00  8.2184E-01  1.0069E+00  8.8052E-01  9.3252E-01  9.4479E-01  2.2716E-01  8.5994E-01  4.7634E-01
             2.3431E+00
 PARAMETER:  5.1148E-02  1.0948E-01 -9.6207E-02  1.0684E-01 -2.7244E-02  3.0134E-02  4.3203E-02 -1.3821E+00 -5.0898E-02 -6.4162E-01
             9.5146E-01
 GRADIENT:   2.4764E+01  2.9773E+00 -5.4156E+00  9.5563E+00  1.2062E+01 -6.8455E+00  8.9845E-01  4.0176E-01  1.9883E+00  8.3782E-01
             1.9310E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1746.29115465108        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      325
 NPARAMETR:  9.4909E-01  9.9727E-01  7.5017E-01  1.0063E+00  8.3280E-01  9.3730E-01  9.7410E-01  2.2187E-01  8.4416E-01  4.4780E-01
             2.3222E+00
 PARAMETER:  4.7745E-02  9.7269E-02 -1.8745E-01  1.0629E-01 -8.2961E-02  3.5246E-02  7.3756E-02 -1.4057E+00 -6.9415E-02 -7.0341E-01
             9.4251E-01
 GRADIENT:  -5.4737E+01 -4.2536E+00 -5.9140E+00 -4.5173E+00  1.0322E+01 -1.1906E+01 -7.6029E-02  3.8246E-01  7.8574E-01  5.5091E-01
             5.6711E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1747.91506558433        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      504
 NPARAMETR:  9.7277E-01  8.0186E-01  6.4679E-01  1.1137E+00  6.8031E-01  9.6166E-01  1.1810E+00  1.5146E-01  7.7520E-01  4.2015E-01
             2.2990E+00
 PARAMETER:  7.2394E-02 -1.2082E-01 -3.3574E-01  2.0770E-01 -2.8520E-01  6.0909E-02  2.6634E-01 -1.7875E+00 -1.5463E-01 -7.6714E-01
             9.3247E-01
 GRADIENT:   5.0493E+00  7.4620E+00 -4.6481E+00  1.5281E+01  2.5366E+00 -1.3866E+00  7.5603E-01  2.3295E-01  1.5048E+00  4.0992E-01
             7.4598E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1748.77445652099        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      681
 NPARAMETR:  9.7877E-01  5.7165E-01  7.8394E-01  1.2587E+00  6.7186E-01  9.4942E-01  1.4702E+00  4.6668E-02  7.2478E-01  5.0937E-01
             2.3207E+00
 PARAMETER:  7.8542E-02 -4.5923E-01 -1.4342E-01  3.3004E-01 -2.9771E-01  4.8093E-02  4.8540E-01 -2.9647E+00 -2.2189E-01 -5.7457E-01
             9.4187E-01
 GRADIENT:   2.9478E+01  1.1154E+01  1.4790E+01  7.5082E+00 -2.6575E+01 -4.3807E+00  1.3505E-01  2.4863E-02  1.6137E-01 -2.8960E-02
             1.3196E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1749.90096707774        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      857
 NPARAMETR:  9.6144E-01  4.3412E-01  8.6683E-01  1.3445E+00  6.8164E-01  9.6449E-01  1.7108E+00  1.3449E-02  7.0321E-01  6.0787E-01
             2.2644E+00
 PARAMETER:  6.0673E-02 -7.3443E-01 -4.2915E-02  3.9603E-01 -2.8326E-01  6.3849E-02  6.3696E-01 -4.2089E+00 -2.5210E-01 -3.9780E-01
             9.1730E-01
 GRADIENT:  -3.7437E+00  5.4854E+00  5.0960E+00  8.4623E+00 -9.9937E+00  2.4617E+00 -3.7500E-01  2.2853E-03 -1.0184E+00  1.3378E+00
            -4.3713E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1750.68689464511        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1032
 NPARAMETR:  9.5728E-01  2.3438E-01  8.9190E-01  1.4568E+00  6.4909E-01  9.4924E-01  2.2545E+00  1.0000E-02  6.7242E-01  6.0480E-01
             2.2714E+00
 PARAMETER:  5.6340E-02 -1.3508E+00 -1.4406E-02  4.7622E-01 -3.3218E-01  4.7904E-02  9.1293E-01 -7.2171E+00 -2.9687E-01 -4.0286E-01
             9.2040E-01
 GRADIENT:  -2.9961E+00  2.7433E-01 -2.1000E+00  1.2847E+01  9.5030E-01 -1.9141E+00 -3.8313E+00  0.0000E+00 -2.6848E+00 -1.5166E+00
            -6.2917E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1751.33073288557        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1208
 NPARAMETR:  9.5783E-01  1.3810E-01  9.3965E-01  1.5189E+00  6.4832E-01  9.5185E-01  3.1566E+00  1.0000E-02  6.6786E-01  6.4873E-01
             2.2874E+00
 PARAMETER:  5.6911E-02 -1.8798E+00  3.7749E-02  5.1797E-01 -3.3337E-01  5.0655E-02  1.2495E+00 -1.0042E+01 -3.0367E-01 -3.3274E-01
             9.2744E-01
 GRADIENT:   3.8717E+00  2.4840E+00  3.4730E+00  1.7709E+01 -7.1678E+00 -6.1814E-02  5.6746E-01  0.0000E+00  1.7382E+00  4.1637E-01
             3.2590E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1753.71259235446        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1384
 NPARAMETR:  9.4954E-01  1.3082E-02  9.0681E-01  1.5721E+00  6.0940E-01  9.5302E-01  9.1228E+00  1.0000E-02  6.4530E-01  6.3736E-01
             2.2681E+00
 PARAMETER:  4.8227E-02 -4.2365E+00  2.1744E-03  5.5243E-01 -3.9528E-01  5.1881E-02  2.3108E+00 -2.3054E+01 -3.3804E-01 -3.5042E-01
             9.1895E-01
 GRADIENT:  -7.4359E+00 -9.0784E-01  3.9965E+00  1.6014E+01 -7.9040E+00  1.2717E+00 -1.5866E+00  0.0000E+00  2.9085E-01  3.3387E-01
            -2.5040E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1753.86323691201        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1574
 NPARAMETR:  9.5132E-01  1.3833E-02  9.0453E-01  1.5719E+00  6.1101E-01  9.5010E-01  9.3069E+00  1.0000E-02  6.4506E-01  6.3405E-01
             2.2709E+00
 PARAMETER:  5.0096E-02 -4.1807E+00 -3.4133E-04  5.5228E-01 -3.9265E-01  4.8817E-02  2.3308E+00 -2.3054E+01 -3.3842E-01 -3.5563E-01
             9.2019E-01
 GRADIENT:  -2.9283E+00 -5.8752E-01 -1.4967E+00  1.5749E+01  1.9529E-02  1.6105E-01 -6.9337E-01  0.0000E+00  6.7343E-02  1.4896E-01
            -2.2317E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1753.89020924600        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1749
 NPARAMETR:  9.5283E-01  1.3905E-02  9.1772E-01  1.5704E+00  6.1689E-01  9.4975E-01  9.2983E+00  1.0000E-02  6.4362E-01  6.3037E-01
             2.2787E+00
 PARAMETER:  5.1685E-02 -4.1755E+00  1.4136E-02  5.5134E-01 -3.8307E-01  4.8442E-02  2.3298E+00 -2.3054E+01 -3.4065E-01 -3.6145E-01
             9.2362E-01
 GRADIENT:   8.2884E-01 -4.5271E-01  9.8668E-02  3.5045E+00  1.9360E-01  1.5351E-01 -3.3249E-01  0.0000E+00 -7.6456E-02 -9.7900E-02
            -1.1644E-01

0ITERATION NO.:   61    OBJECTIVE VALUE:  -1753.89020924600        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1771
 NPARAMETR:  9.5283E-01  1.3905E-02  9.1772E-01  1.5704E+00  6.1689E-01  9.4975E-01  9.2983E+00  1.0000E-02  6.4362E-01  6.3037E-01
             2.2787E+00
 PARAMETER:  5.1685E-02 -4.1755E+00  1.4136E-02  5.5134E-01 -3.8307E-01  4.8442E-02  2.3298E+00 -2.3054E+01 -3.4065E-01 -3.6145E-01
             9.2362E-01
 GRADIENT:   8.2884E-01 -4.5271E-01  9.8668E-02  3.5045E+00  1.9360E-01  1.5351E-01 -3.3249E-01  0.0000E+00 -7.6456E-02 -9.7900E-02
            -1.1644E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1771
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.5701E-04  7.0859E-03  6.3155E-05 -1.0426E-02 -1.2226E-02
 SE:             2.9377E-02  5.6428E-03  1.8278E-04  2.7408E-02  1.7832E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8759E-01  2.0921E-01  7.2970E-01  7.0365E-01  4.9296E-01

 ETASHRINKSD(%)  1.5819E+00  8.1096E+01  9.9388E+01  8.1803E+00  4.0260E+01
 ETASHRINKVR(%)  3.1388E+00  9.6426E+01  9.9996E+01  1.5691E+01  6.4311E+01
 EBVSHRINKSD(%)  1.6449E+00  8.7118E+01  9.9346E+01  8.1876E+00  4.0827E+01
 EBVSHRINKVR(%)  3.2627E+00  9.8341E+01  9.9996E+01  1.5705E+01  6.4986E+01
 RELATIVEINF(%)  9.6328E+01  1.0053E+00  1.9985E-04  3.0783E+01  1.6457E+00
 EPSSHRINKSD(%)  2.4056E+01
 EPSSHRINKVR(%)  4.2325E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1753.8902092459966     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -834.95167604132394     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.63
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.41
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1753.890       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.53E-01  1.39E-02  9.18E-01  1.57E+00  6.17E-01  9.50E-01  9.30E+00  1.00E-02  6.44E-01  6.30E-01  2.28E+00
 


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
+        1.31E+03
 
 TH 2
+       -5.72E+02  1.03E+06
 
 TH 3
+       -1.91E+01  2.04E+03  5.52E+02
 
 TH 4
+       -1.66E+01 -3.69E+04 -1.86E+02  2.70E+03
 
 TH 5
+        4.89E+01  1.47E+05 -1.12E+03  6.41E+01  2.52E+03
 
 TH 6
+        2.50E+00 -4.38E+01  2.62E+00 -9.31E+00 -2.09E+00  2.05E+02
 
 TH 7
+       -1.37E+00  2.71E+03  5.21E+00 -1.00E+02 -1.45E+01 -9.72E-02  1.04E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.03E+00  2.94E+02  2.68E+01 -6.48E+01 -2.20E+01 -1.57E+00  8.53E-01  0.00E+00  4.13E+02
 
 TH10
+       -5.39E+00 -1.56E+05 -8.69E+00 -7.12E+01 -1.18E+01  5.58E-01  5.26E+00  0.00E+00  1.34E+01  6.43E+01
 
 TH11
+       -6.09E+00 -1.62E+04 -3.59E+01  7.42E+02  7.75E+01  3.73E+00 -4.35E+01  0.00E+00 -1.95E+00  1.45E+00  4.16E+02
 
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
 #CPUT: Total CPU Time in Seconds,       33.112
Stop Time:
Wed Sep 29 23:51:05 CDT 2021
