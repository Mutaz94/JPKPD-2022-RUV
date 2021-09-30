Wed Sep 29 13:35:11 CDT 2021
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
$DATA ../../../../data/spa/A3/dat49.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   25.1395980985696        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.6900E+02  7.5605E+01  9.6306E+01  5.9316E+00  1.6451E+02  3.5601E+01 -1.3064E+02 -4.1323E+01 -2.0448E+02 -1.6876E+02
            -2.7032E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -960.280526772306        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0801E+00  8.9455E-01  7.8503E-01  1.3101E+00  7.6400E-01  7.6928E-01  1.6020E+00  1.0092E+00  1.9517E+00  1.5177E+00
             8.8251E+00
 PARAMETER:  1.7709E-01 -1.1435E-02 -1.4203E-01  3.7013E-01 -1.6919E-01 -1.6230E-01  5.7123E-01  1.0913E-01  7.6873E-01  5.1718E-01
             2.2776E+00
 GRADIENT:   3.7084E+01 -1.9454E+01 -2.2856E+01 -1.9025E+00 -1.1292E+01 -8.5830E+00  1.4545E+01  6.4450E+00  6.5334E+01  3.7865E+01
             3.5108E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1015.73194459758        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.7841E-01  2.8817E-01  5.7708E-01  1.9105E+00  3.6090E-01  9.0572E-01  1.0535E+01  5.3852E-01  1.4979E+00  7.7496E-01
             6.4449E+00
 PARAMETER:  7.8174E-02 -1.1442E+00 -4.4978E-01  7.4734E-01 -9.1915E-01  9.7934E-04  2.4547E+00 -5.1894E-01  5.0406E-01 -1.5495E-01
             1.9633E+00
 GRADIENT:  -7.0535E+01  3.0141E+01  3.0453E+01  1.3142E+02 -8.4640E+01  6.7613E-01  3.4275E+01  2.1378E+00  3.1740E+01  4.0136E+00
             2.5460E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1171.12783581248        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  8.8755E-01  4.6523E-01  1.8478E-01  9.9989E-01  2.4686E-01  1.0542E+00  1.0927E+00  5.9546E-02  1.3208E+00  5.6682E-01
             3.4926E+00
 PARAMETER: -1.9285E-02 -6.6523E-01 -1.5886E+00  9.9889E-02 -1.2989E+00  1.5279E-01  1.8862E-01 -2.7210E+00  3.7825E-01 -4.6772E-01
             1.3506E+00
 GRADIENT:  -1.1363E+02  4.9174E+01  9.7668E+01 -6.8098E+01 -7.7851E+01  2.0232E+01 -8.4845E+00  1.2783E-03 -3.1876E+01 -1.3203E+01
             5.5138E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1187.54584229973        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      386
 NPARAMETR:  9.6214E-01  4.5682E-01  1.9100E-01  1.0845E+00  2.5704E-01  9.7457E-01  1.1069E+00  1.6387E-02  1.4745E+00  6.4051E-01
             3.0568E+00
 PARAMETER:  6.1405E-02 -6.8346E-01 -1.5555E+00  1.8113E-01 -1.2585E+00  7.4238E-02  2.0160E-01 -4.0113E+00  4.8834E-01 -3.4548E-01
             1.2174E+00
 GRADIENT:   3.3256E+01  1.3407E+01  4.0828E+01 -5.4396E+00 -4.3764E+01 -3.2254E+00 -9.1459E+00  4.2365E-04 -6.8732E+00 -1.0026E+01
            -2.2145E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1193.82166855489        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      562
 NPARAMETR:  9.2914E-01  4.7090E-01  1.3223E-01  1.0134E+00  2.2700E-01  1.0237E+00  1.0044E+00  1.0000E-02  1.6203E+00  7.0587E-01
             2.8928E+00
 PARAMETER:  2.6509E-02 -6.5312E-01 -1.9232E+00  1.1331E-01 -1.3828E+00  1.2347E-01  1.0439E-01 -5.5112E+00  5.8259E-01 -2.4832E-01
             1.1622E+00
 GRADIENT:  -1.9027E+00  1.0710E+01  2.8223E+00  5.2415E+00 -1.2582E+01  7.4384E+00  5.1910E-01  0.0000E+00 -2.2278E+01 -5.9048E+00
            -1.5084E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1194.38039054704        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      682
 NPARAMETR:  9.2103E-01  4.5695E-01  1.3273E-01  9.9614E-01  2.2827E-01  9.8579E-01  8.6517E-01  1.0000E-02  1.6307E+00  7.8747E-01
             2.9067E+00
 PARAMETER:  1.7740E-02 -6.8319E-01 -1.9194E+00  9.6134E-02 -1.3772E+00  8.5688E-02 -4.4831E-02 -5.5600E+00  5.8901E-01 -1.3893E-01
             1.1670E+00
 GRADIENT:   9.7202E+00  3.9577E-01  2.6174E+01  4.6519E-01  1.3605E+02 -3.9111E-01  1.6721E+00  0.0000E+00 -2.0740E+00  4.1518E-01
            -7.9393E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1194.53216402556        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      753
 NPARAMETR:  9.1451E-01  4.5350E-01  1.3255E-01  9.8977E-01  2.2738E-01  9.8434E-01  4.8051E-01  1.0000E-02  1.6306E+00  8.6842E-01
             2.9074E+00
 PARAMETER:  1.0629E-02 -6.9076E-01 -1.9208E+00  8.9716E-02 -1.3811E+00  8.4216E-02 -6.3291E-01 -5.5600E+00  5.8892E-01 -4.1078E-02
             1.1673E+00
 GRADIENT:  -8.2999E+00  1.4582E+00  3.4300E+01 -7.1054E-02  1.3693E+02  2.6608E-01  5.2439E-01  0.0000E+00 -4.9294E+00  1.8087E+00
            -1.6415E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1196.15540427156        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      827
 NPARAMETR:  9.3228E-01  3.8921E-01  1.2732E-01  9.8791E-01  2.0453E-01  9.5443E-01  1.0000E-02  1.0000E-02  1.6275E+00  9.1919E-01
             2.9249E+00
 PARAMETER:  2.9880E-02 -8.4364E-01 -1.9611E+00  8.7840E-02 -1.4870E+00  5.3360E-02 -1.2927E+01 -5.5600E+00  5.8704E-01  1.5734E-02
             1.1733E+00
 GRADIENT:   1.6842E+01  7.4383E+00  5.0240E+01  4.2368E+00  1.1685E+02 -8.8917E+00  0.0000E+00  0.0000E+00 -5.7930E+00  6.1565E+00
            -1.3535E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1196.52793451152        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1006
 NPARAMETR:  9.3775E-01  3.9440E-01  1.2743E-01  9.9336E-01  2.0522E-01  9.9116E-01  1.0000E-02  1.0000E-02  1.6277E+00  8.9375E-01
             2.9248E+00
 PARAMETER:  3.5725E-02 -8.3039E-01 -1.9602E+00  9.3340E-02 -1.4837E+00  9.1122E-02 -1.2562E+01 -5.5600E+00  5.8715E-01 -1.2324E-02
             1.1732E+00
 GRADIENT:   1.0947E-01  5.9384E-02  2.0646E+01  2.6024E-01 -8.5475E+00  6.3896E-01  0.0000E+00  0.0000E+00 -2.3981E+01 -2.3085E-01
            -1.2435E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1197.71917554216        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1187
 NPARAMETR:  9.3594E-01  4.0148E-01  1.1666E-01  9.7793E-01  2.0148E-01  9.9441E-01  1.0000E-02  1.0000E-02  1.6405E+00  8.9593E-01
             3.0246E+00
 PARAMETER:  3.3794E-02 -8.1259E-01 -2.0485E+00  7.7681E-02 -1.5021E+00  9.4394E-02 -1.6424E+01 -5.5600E+00  5.9497E-01 -9.8959E-03
             1.2068E+00
 GRADIENT:   2.8090E+00 -8.5287E-01  1.2205E+00  3.1855E+00  5.9582E+00  6.8305E-01  0.0000E+00  0.0000E+00 -3.1513E+01  5.9852E-01
             1.1688E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1198.25374387121        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1362
 NPARAMETR:  9.3501E-01  3.9363E-01  1.1548E-01  9.7042E-01  1.9860E-01  9.9226E-01  1.0000E-02  1.0000E-02  1.6634E+00  8.9484E-01
             3.0364E+00
 PARAMETER:  3.2945E-02 -8.3254E-01 -2.0632E+00  6.8970E-02 -1.5199E+00  9.2178E-02 -1.8641E+01 -5.5600E+00  6.1023E-01 -1.1191E-02
             1.2080E+00
 GRADIENT:   6.7707E-01 -2.2736E-01 -1.5284E+03 -7.1892E-01 -2.0956E+03 -4.2820E-02  0.0000E+00  0.0000E+00  5.1881E+03 -4.0132E-02
            -2.6851E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1362
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7313E-03 -5.8043E-05  2.8565E-04 -1.7613E-02  6.9055E-03
 SE:             2.8559E-02  1.3843E-04  1.3382E-04  2.7664E-02  2.5866E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5166E-01  6.7501E-01  3.2798E-02  5.2435E-01  7.8949E-01

 ETASHRINKSD(%)  4.3239E+00  9.9536E+01  9.9552E+01  7.3224E+00  1.3345E+01
 ETASHRINKVR(%)  8.4609E+00  9.9998E+01  9.9998E+01  1.4109E+01  2.4909E+01
 EBVSHRINKSD(%)  3.7562E+00  9.9509E+01  9.9553E+01  1.0709E+01  1.4117E+01
 EBVSHRINKVR(%)  7.3713E+00  9.9998E+01  9.9998E+01  2.0271E+01  2.6242E+01
 RELATIVEINF(%)  8.7323E+01  2.5637E-04  3.5974E-04  5.2806E+01  3.9539E+00
 EPSSHRINKSD(%)  3.5449E+01
 EPSSHRINKVR(%)  5.8332E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1198.2537438712143     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -463.10291730747610     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.82
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.38
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1198.254       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.35E-01  3.94E-01  1.15E-01  9.69E-01  1.98E-01  9.92E-01  1.00E-02  1.00E-02  1.67E+00  8.95E-01  3.03E+00
 


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
+        1.20E+03
 
 TH 2
+        6.07E+01  1.76E+03
 
 TH 3
+       -1.31E+03  3.73E+03  1.41E+06
 
 TH 4
+        2.54E+00  2.96E+01  1.59E+03  3.43E+02
 
 TH 5
+        2.83E+06 -8.16E+05 -1.94E+04  1.61E+03  9.05E+05
 
 TH 6
+        3.09E+00 -1.54E+01 -8.83E+02 -1.14E+01 -7.56E+02  1.73E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.77E+02  2.39E+05  1.40E+03  8.07E+05  8.48E+02 -7.89E+05  0.00E+00  0.00E+00  7.69E+04
 
 TH10
+       -4.04E+00 -1.97E+01 -9.21E+02  8.15E+00 -6.88E+02  6.02E+00  0.00E+00  0.00E+00  2.35E+02  1.36E+02
 
 TH11
+       -6.75E+01  6.80E+04  3.96E+03  1.46E+02  1.45E+05 -6.36E+01  0.00E+00  0.00E+00 -2.19E+04 -5.53E+01  6.25E+03
 
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
 #CPUT: Total CPU Time in Seconds,       25.258
Stop Time:
Wed Sep 29 13:35:38 CDT 2021
