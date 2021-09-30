Wed Sep 29 11:16:38 CDT 2021
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
$DATA ../../../../data/spa/B/dat47.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m47.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1643.55506350210        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3390E+02 -3.8599E+01 -7.0979E+00 -4.9003E+01  3.6467E+01  6.1981E+01 -3.7716E+00 -1.6730E-01 -2.4344E+01 -2.1105E+00
             9.4785E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1650.56995186089        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.8462E-01  1.0134E+00  1.0100E+00  1.0597E+00  9.7529E-01  9.2931E-01  1.0145E+00  1.0056E+00  1.0918E+00  1.0059E+00
             9.7641E-01
 PARAMETER:  8.4501E-02  1.1334E-01  1.0995E-01  1.5802E-01  7.4984E-02  2.6692E-02  1.1439E-01  1.0557E-01  1.8786E-01  1.0587E-01
             7.6126E-02
 GRADIENT:   3.7192E+00 -1.4085E+00  1.1747E+00 -9.2180E+00 -5.2490E+00 -1.5727E+00  1.1192E+00 -7.6268E-01 -2.2833E+00  4.2898E-01
            -2.4074E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1650.82374882828        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      353
 NPARAMETR:  9.8297E-01  8.8314E-01  1.2702E+00  1.1685E+00  1.0298E+00  9.2506E-01  8.2547E-01  1.1628E+00  1.1191E+00  1.0781E+00
             9.7253E-01
 PARAMETER:  8.2819E-02 -2.4266E-02  3.3917E-01  2.5572E-01  1.2934E-01  2.2103E-02 -9.1798E-02  2.5082E-01  2.1256E-01  1.7520E-01
             7.2146E-02
 GRADIENT:   5.1145E+00  1.3848E+01  1.2686E+01  1.7256E+01 -1.3666E+01 -2.8467E+00  4.6210E-01 -4.8689E+00  8.4574E+00 -4.3945E+00
            -4.6343E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1651.32961168573        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.8142E-01  8.9450E-01  1.2749E+00  1.1537E+00  1.0473E+00  9.3165E-01  8.2022E-01  1.2364E+00  1.0944E+00  1.1126E+00
             9.7985E-01
 PARAMETER:  8.1241E-02 -1.1492E-02  3.4290E-01  2.4297E-01  1.4624E-01  2.9197E-02 -9.8181E-02  3.1224E-01  1.9021E-01  2.0672E-01
             7.9648E-02
 GRADIENT:   8.5516E-01  6.4002E+00  2.2341E+00  6.3059E+00 -4.8769E+00  2.9646E-02  2.6309E-01 -1.4552E-01  4.4505E-01  2.9803E-01
            -2.3614E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1651.61044926423        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  9.7832E-01  6.3943E-01  1.3892E+00  1.3154E+00  1.0013E+00  9.3143E-01  7.9369E-01  1.2405E+00  9.8700E-01  1.1012E+00
             9.8068E-01
 PARAMETER:  7.8081E-02 -3.4718E-01  4.2872E-01  3.7414E-01  1.0130E-01  2.8969E-02 -1.3107E-01  3.1548E-01  8.6911E-02  1.9641E-01
             8.0496E-02
 GRADIENT:  -5.3134E-01  2.7799E+00  1.7614E+00  2.7032E+00 -1.9601E+00  4.9919E-01 -4.5106E-01 -7.3940E-01 -8.3848E-01 -6.8350E-02
             1.4845E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1651.61444546167        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  9.7784E-01  5.8494E-01  1.4055E+00  1.3531E+00  9.8826E-01  9.3011E-01  8.0566E-01  1.2415E+00  9.6531E-01  1.0941E+00
             9.8020E-01
 PARAMETER:  7.7588E-02 -4.3625E-01  4.4037E-01  4.0237E-01  8.8189E-02  2.7548E-02 -1.1609E-01  3.1629E-01  6.4697E-02  1.8992E-01
             7.9999E-02
 GRADIENT:  -2.0589E-01  3.8691E+00  2.3099E+00  7.0502E+00 -2.7041E+00  4.8843E-02 -4.8247E-01 -1.0561E+00 -7.4595E-01 -1.7636E-01
            -7.3540E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1651.61802151953        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  9.7726E-01  5.2486E-01  1.4201E+00  1.3933E+00  9.7360E-01  9.2875E-01  8.4466E-01  1.2484E+00  9.4126E-01  1.0842E+00
             9.7947E-01
 PARAMETER:  7.6999E-02 -5.4462E-01  4.5070E-01  4.3169E-01  7.3243E-02  2.6090E-02 -6.8821E-02  3.2185E-01  3.9463E-02  1.8089E-01
             7.9259E-02
 GRADIENT:   1.6751E-01  4.1731E+00  2.3110E+00  1.0008E+01 -2.9069E+00 -4.0325E-01 -4.3265E-01 -1.1050E+00 -5.9971E-01 -2.5408E-01
            -2.9819E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1651.61952725214        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1229
 NPARAMETR:  9.7659E-01  4.6395E-01  1.4349E+00  1.4332E+00  9.5956E-01  9.2761E-01  9.1824E-01  1.2621E+00  9.1701E-01  1.0738E+00
             9.7876E-01
 PARAMETER:  7.6314E-02 -6.6798E-01  4.6107E-01  4.5991E-01  5.8718E-02  2.4853E-02  1.4708E-02  3.3281E-01  1.3364E-02  1.7122E-01
             7.8532E-02
 GRADIENT:   4.8326E-01  3.9380E+00  1.9811E+00  1.1358E+01 -2.7283E+00 -7.6753E-01 -3.3021E-01 -9.7841E-01 -4.1018E-01 -2.8424E-01
            -4.6793E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1651.62174115035        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1404
 NPARAMETR:  9.7606E-01  4.2456E-01  1.4467E+00  1.4582E+00  9.5176E-01  9.2731E-01  9.8668E-01  1.2759E+00  9.0156E-01  1.0680E+00
             9.7850E-01
 PARAMETER:  7.5769E-02 -7.5670E-01  4.6929E-01  4.7717E-01  5.0552E-02  2.4538E-02  8.6585E-02  3.4369E-01 -3.6340E-03  1.6577E-01
             7.8265E-02
 GRADIENT:   5.2514E-01  3.4695E+00  1.6684E+00  1.0613E+01 -2.4154E+00 -8.1149E-01 -2.4476E-01 -8.2179E-01 -2.7641E-01 -2.4796E-01
            -4.7382E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1651.62200081721        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1584
 NPARAMETR:  9.7593E-01  4.1615E-01  1.4498E+00  1.4634E+00  9.5033E-01  9.2732E-01  1.0019E+00  1.2795E+00  8.9833E-01  1.0670E+00
             9.7848E-01
 PARAMETER:  7.5638E-02 -7.7672E-01  4.7140E-01  4.8075E-01  4.9052E-02  2.4548E-02  1.0193E-01  3.4650E-01 -7.2215E-03  1.6483E-01
             7.8248E-02
 GRADIENT:   5.0882E-01  3.3402E+00  1.6019E+00  1.0251E+01 -2.3291E+00 -7.8972E-01 -2.2796E-01 -7.8678E-01 -2.5277E-01 -2.3496E-01
            -4.5901E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1651.69224783762        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1772             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7614E-01  4.0794E-01  1.4495E+00  1.4587E+00  9.4975E-01  9.2934E-01  1.0708E+00  1.2878E+00  8.9275E-01  1.0663E+00
             9.7872E-01
 PARAMETER:  7.5853E-02 -7.9664E-01  4.7124E-01  4.7756E-01  4.8443E-02  2.6724E-02  1.6838E-01  3.5291E-01 -1.3449E-02  1.6417E-01
             7.8491E-02
 GRADIENT:   4.1116E+02  5.2219E+01  9.5413E+00  6.0751E+02  7.2320E+00  3.9491E+01  9.8262E-01  1.6395E+00  1.4233E+01  1.5185E+00
             7.3648E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1651.69312537957        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1958             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7614E-01  4.0836E-01  1.4486E+00  1.4588E+00  9.4927E-01  9.2934E-01  1.0761E+00  1.2876E+00  8.9358E-01  1.0656E+00
             9.7875E-01
 PARAMETER:  7.5854E-02 -7.9562E-01  4.7061E-01  4.7759E-01  4.7937E-02  2.6721E-02  1.7332E-01  3.5281E-01 -1.2517E-02  1.6356E-01
             7.8519E-02
 GRADIENT:   4.1110E+02  5.2329E+01  9.6224E+00  6.0795E+02  6.9648E+00  3.9483E+01  1.0455E+00  1.6804E+00  1.4573E+01  1.5278E+00
             7.6605E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1651.69372448741        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2143            RESET HESSIAN, TYPE II
 NPARAMETR:  9.7615E-01  4.0855E-01  1.4474E+00  1.4584E+00  9.4883E-01  9.2934E-01  1.0838E+00  1.2865E+00  8.9367E-01  1.0650E+00
             9.7873E-01
 PARAMETER:  7.5865E-02 -7.9513E-01  4.6976E-01  4.7735E-01  4.7470E-02  2.6725E-02  1.8049E-01  3.5193E-01 -1.2419E-02  1.6297E-01
             7.8498E-02
 GRADIENT:   4.1112E+02  5.2260E+01  9.6439E+00  6.0733E+02  6.9607E+00  3.9481E+01  1.1245E+00  1.6699E+00  1.4663E+01  1.5116E+00
             7.6853E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1651.69417058765        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2326
 NPARAMETR:  9.7616E-01  4.0921E-01  1.4461E+00  1.4581E+00  9.4845E-01  9.2934E-01  1.0815E+00  1.2855E+00  8.9393E-01  1.0646E+00
             9.7873E-01
 PARAMETER:  7.5873E-02 -7.9354E-01  4.6888E-01  4.7714E-01  4.7071E-02  2.6723E-02  1.7839E-01  3.5112E-01 -1.2128E-02  1.6261E-01
             7.8501E-02
 GRADIENT:   1.5043E+00  1.0076E-01  2.9690E-01 -6.9292E+00 -3.2809E-02  1.0385E-01 -1.3461E-02  6.8844E-03 -1.0930E-01  1.3615E-02
            -3.0098E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1651.69478051721        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2511             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7618E-01  4.0966E-01  1.4443E+00  1.4575E+00  9.4817E-01  9.2935E-01  1.0882E+00  1.2845E+00  8.9453E-01  1.0643E+00
             9.7876E-01
 PARAMETER:  7.5889E-02 -7.9243E-01  4.6764E-01  4.7673E-01  4.6783E-02  2.6726E-02  1.8452E-01  3.5033E-01 -1.1452E-02  1.6227E-01
             7.8536E-02
 GRADIENT:   4.1106E+02  5.2250E+01  9.4495E+00  6.0564E+02  7.1274E+00  3.9470E+01  1.1839E+00  1.7111E+00  1.4881E+01  1.5195E+00
             8.1579E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1651.69507603180        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2696
 NPARAMETR:  9.7619E-01  4.1010E-01  1.4431E+00  1.4571E+00  9.4795E-01  9.2935E-01  1.0922E+00  1.2835E+00  8.9469E-01  1.0640E+00
             9.7877E-01
 PARAMETER:  7.5901E-02 -7.9134E-01  4.6678E-01  4.7642E-01  4.6545E-02  2.6729E-02  1.8824E-01  3.4963E-01 -1.1279E-02  1.6201E-01
             7.8546E-02
 GRADIENT:   1.5391E+00 -1.0606E-01  2.1483E-03 -7.6336E+00  3.5999E-01  1.0364E-01  2.5560E-02  7.8377E-02  1.8181E-01  3.4512E-02
             3.3307E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1651.69564734754        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2881             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7620E-01  4.1108E-01  1.4420E+00  1.4567E+00  9.4753E-01  9.2935E-01  1.0818E+00  1.2821E+00  8.9467E-01  1.0636E+00
             9.7872E-01
 PARAMETER:  7.5909E-02 -7.8896E-01  4.6604E-01  4.7619E-01  4.6100E-02  2.6733E-02  1.7864E-01  3.4847E-01 -1.1306E-02  1.6165E-01
             7.8487E-02
 GRADIENT:   4.1106E+02  5.2491E+01  9.5465E+00  6.0437E+02  6.8792E+00  3.9474E+01  1.1195E+00  1.6604E+00  1.4665E+01  1.4991E+00
             7.8154E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1651.69585574790        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3063
 NPARAMETR:  9.7620E-01  4.1174E-01  1.4412E+00  1.4564E+00  9.4726E-01  9.2935E-01  1.0789E+00  1.2810E+00  8.9476E-01  1.0633E+00
             9.7869E-01
 PARAMETER:  7.5916E-02 -7.8737E-01  4.6544E-01  4.7597E-01  4.5814E-02  2.6735E-02  1.7597E-01  3.4767E-01 -1.1197E-02  1.6141E-01
             7.8464E-02
 GRADIENT:   1.4935E+00  1.4996E-01  2.8267E-01 -6.6817E+00 -1.8053E-01  1.0280E-01 -2.1082E-02 -4.4805E-03 -1.5831E-01  9.9157E-03
            -3.6461E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1651.69632315827        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     3252             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7623E-01  4.1215E-01  1.4390E+00  1.4555E+00  9.4716E-01  9.2936E-01  1.0919E+00  1.2801E+00  8.9549E-01  1.0630E+00
             9.7876E-01
 PARAMETER:  7.5942E-02 -7.8637E-01  4.6397E-01  4.7534E-01  4.5711E-02  2.6742E-02  1.8793E-01  3.4694E-01 -1.0381E-02  1.6107E-01
             7.8533E-02
 GRADIENT:   4.1104E+02  5.2309E+01  9.2123E+00  6.0189E+02  7.3624E+00  3.9463E+01  1.2419E+00  1.7144E+00  1.4932E+01  1.4988E+00
             8.4298E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1651.69660961796        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3439             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7623E-01  4.1253E-01  1.4383E+00  1.4554E+00  9.4693E-01  9.2936E-01  1.0917E+00  1.2795E+00  8.9555E-01  1.0628E+00
             9.7876E-01
 PARAMETER:  7.5944E-02 -7.8544E-01  4.6346E-01  4.7527E-01  4.5467E-02  2.6741E-02  1.8771E-01  3.4646E-01 -1.0320E-02  1.6094E-01
             7.8527E-02
 GRADIENT:   4.1102E+02  5.2399E+01  9.2259E+00  6.0185E+02  7.2664E+00  3.9462E+01  1.2397E+00  1.7134E+00  1.4903E+01  1.5129E+00
             8.4119E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1651.69700376926        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3626
 NPARAMETR:  9.7624E-01  4.1314E-01  1.4375E+00  1.4549E+00  9.4671E-01  9.2937E-01  1.0882E+00  1.2784E+00  8.9575E-01  1.0626E+00
             9.7874E-01
 PARAMETER:  7.5958E-02 -7.8397E-01  4.6287E-01  4.7492E-01  4.5237E-02  2.6747E-02  1.8454E-01  3.4563E-01 -1.0095E-02  1.6068E-01
             7.8507E-02
 GRADIENT:   1.5454E+00 -1.1034E-01 -4.1572E-02 -7.6056E+00  2.8477E-01  1.0497E-01  1.6748E-02  5.9568E-02  1.3588E-01  1.5342E-02
             2.6845E-02

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1651.69738867584        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3811             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7625E-01  4.1408E-01  1.4368E+00  1.4546E+00  9.4645E-01  9.2937E-01  1.0822E+00  1.2774E+00  8.9575E-01  1.0623E+00
             9.7869E-01
 PARAMETER:  7.5963E-02 -7.8170E-01  4.6240E-01  4.7472E-01  4.4964E-02  2.6750E-02  1.7898E-01  3.4479E-01 -1.0089E-02  1.6048E-01
             7.8459E-02
 GRADIENT:   4.1105E+02  5.2698E+01  9.4672E+00  6.0058E+02  6.8505E+00  3.9470E+01  1.1434E+00  1.6279E+00  1.4669E+01  1.4823E+00
             7.8581E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1651.69755915179        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3996
 NPARAMETR:  9.7626E-01  4.1468E-01  1.4362E+00  1.4543E+00  9.4630E-01  9.2937E-01  1.0803E+00  1.2766E+00  8.9589E-01  1.0622E+00
             9.7868E-01
 PARAMETER:  7.5971E-02 -7.8025E-01  4.6199E-01  4.7452E-01  4.4809E-02  2.6752E-02  1.7727E-01  3.4423E-01 -9.9417E-03  1.6035E-01
             7.8445E-02
 GRADIENT:   1.5029E+00  1.1335E-01  2.0376E-01 -6.7648E+00 -1.8030E-01  1.0308E-01 -1.2090E-02 -2.8520E-03 -9.5943E-02  9.0587E-03
            -2.5688E-02

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1651.69780267087        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     4169
 NPARAMETR:  9.7627E-01  4.1493E-01  1.4356E+00  1.4540E+00  9.4629E-01  9.2937E-01  1.0826E+00  1.2763E+00  8.9612E-01  1.0621E+00
             9.7869E-01
 PARAMETER:  7.5979E-02 -7.7963E-01  4.6157E-01  4.7429E-01  4.4793E-02  2.6755E-02  1.7940E-01  3.4398E-01 -9.6854E-03  1.6028E-01
             7.8462E-02
 GRADIENT:  -1.5940E-02 -3.3092E-02  1.3872E-01  4.0193E-01 -1.8426E-02 -5.4839E-04 -2.5135E-03  1.2041E-02 -6.6188E-02  7.6694E-03
            -6.2800E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     4169
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0215E-04 -1.0676E-02 -3.6566E-02 -2.7943E-03 -4.1520E-02
 SE:             2.9838E-02  7.0073E-03  1.7875E-02  2.8651E-02  2.1016E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9192E-01  1.2761E-01  4.0794E-02  9.2231E-01  4.8193E-02

 ETASHRINKSD(%)  3.8934E-02  7.6525E+01  4.0116E+01  4.0148E+00  2.9595E+01
 ETASHRINKVR(%)  7.7853E-02  9.4489E+01  6.4140E+01  7.8684E+00  5.0432E+01
 EBVSHRINKSD(%)  4.4742E-01  7.7215E+01  4.4098E+01  4.2244E+00  2.6050E+01
 EBVSHRINKVR(%)  8.9284E-01  9.4808E+01  6.8750E+01  8.2704E+00  4.5314E+01
 RELATIVEINF(%)  9.6820E+01  2.3055E-01  8.4510E+00  5.0121E+00  8.2399E+00
 EPSSHRINKSD(%)  4.5709E+01
 EPSSHRINKVR(%)  7.0525E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1651.6978026708732     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -916.54697610713504     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    56.11
 Elapsed covariance  time in seconds:     5.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1651.698       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  4.15E-01  1.44E+00  1.45E+00  9.46E-01  9.29E-01  1.08E+00  1.28E+00  8.96E-01  1.06E+00  9.79E-01
 


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
 
         2.93E-02  4.84E-01  3.34E-01  2.99E-01  1.86E-01  6.55E-02  1.35E+00  2.93E-01  1.94E-01  2.49E-01  6.45E-02
 


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
+        8.61E-04
 
 TH 2
+        5.60E-03  2.35E-01
 
 TH 3
+       -9.09E-04 -3.50E-02  1.11E-01
 
 TH 4
+       -3.37E-03 -1.43E-01  2.62E-02  8.94E-02
 
 TH 5
+        1.13E-03  6.03E-02  3.35E-02 -3.48E-02  3.44E-02
 
 TH 6
+       -7.35E-05 -3.30E-03  3.28E-03  1.93E-03  4.84E-04  4.28E-03
 
 TH 7
+       -1.26E-02 -6.16E-01  7.36E-02  3.71E-01 -1.68E-01  7.61E-03  1.84E+00
 
 TH 8
+       -1.64E-03 -4.69E-02  7.11E-02  3.23E-02  1.26E-02  3.60E-03  1.27E-01  8.58E-02
 
 TH 9
+        2.36E-03  8.63E-02 -1.91E-02 -5.26E-02  1.99E-02 -1.21E-03 -2.16E-01 -2.15E-02  3.77E-02
 
 TH10
+        9.29E-04  6.61E-02  4.32E-02 -3.83E-02  4.06E-02  2.46E-03 -1.83E-01  1.16E-02  2.11E-02  6.21E-02
 
 TH11
+        1.12E-04  8.35E-04  3.03E-03 -1.92E-04  1.43E-03  1.84E-04 -4.84E-03  3.34E-04 -2.34E-04  1.36E-03  4.16E-03
 
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
+        2.93E-02
 
 TH 2
+        3.94E-01  4.84E-01
 
 TH 3
+       -9.28E-02 -2.16E-01  3.34E-01
 
 TH 4
+       -3.84E-01 -9.90E-01  2.63E-01  2.99E-01
 
 TH 5
+        2.08E-01  6.71E-01  5.41E-01 -6.27E-01  1.86E-01
 
 TH 6
+       -3.82E-02 -1.04E-01  1.50E-01  9.85E-02  3.99E-02  6.55E-02
 
 TH 7
+       -3.18E-01 -9.38E-01  1.63E-01  9.17E-01 -6.67E-01  8.58E-02  1.35E+00
 
 TH 8
+       -1.91E-01 -3.30E-01  7.27E-01  3.69E-01  2.32E-01  1.87E-01  3.20E-01  2.93E-01
 
 TH 9
+        4.14E-01  9.18E-01 -2.95E-01 -9.06E-01  5.54E-01 -9.53E-02 -8.22E-01 -3.78E-01  1.94E-01
 
 TH10
+        1.27E-01  5.48E-01  5.19E-01 -5.15E-01  8.79E-01  1.51E-01 -5.42E-01  1.59E-01  4.36E-01  2.49E-01
 
 TH11
+        5.94E-02  2.67E-02  1.41E-01 -9.95E-03  1.19E-01  4.37E-02 -5.53E-02  1.77E-02 -1.87E-02  8.43E-02  6.45E-02
 
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
+        1.60E+03
 
 TH 2
+       -2.03E+02  5.27E+02
 
 TH 3
+       -1.24E+02  6.88E+01  1.24E+02
 
 TH 4
+       -1.13E+02  5.44E+02 -9.54E-01  7.39E+02
 
 TH 5
+        1.86E+02 -1.72E+02 -2.21E+02 -4.43E+01  6.72E+02
 
 TH 6
+       -5.62E+01  6.83E+01  2.11E+01  5.72E+01  1.08E+01  2.72E+02
 
 TH 7
+       -1.63E+01  3.05E+01  3.26E-01  1.95E+01  6.77E+00  4.53E+00  5.92E+00
 
 TH 8
+        4.40E+01 -3.46E+01 -2.75E+01 -2.90E+01 -1.25E+00 -2.18E+01 -3.61E+00  3.38E+01
 
 TH 9
+       -5.38E+01 -1.43E+02  2.46E+01 -8.66E+01 -5.60E+01 -2.12E+01 -1.12E+01  1.04E+01  2.13E+02
 
 TH10
+        5.19E+01 -1.38E+01 -1.62E+01 -4.83E-02 -9.38E+01 -4.40E+01 -3.07E+00  1.82E+01  1.22E+01  8.75E+01
 
 TH11
+       -2.07E+01 -3.46E+01 -1.67E+01 -3.85E+01 -6.78E+00 -2.05E+01  2.68E-03  1.35E+01  2.26E+01  1.44E+01  2.57E+02
 
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
 #CPUT: Total CPU Time in Seconds,       61.882
Stop Time:
Wed Sep 29 11:17:50 CDT 2021
