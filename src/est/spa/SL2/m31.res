Sat Sep 18 12:11:52 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat31.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m31.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1640.85452538121        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0257E+02 -4.5151E+01 -4.9440E+01 -1.4003E+00  9.0613E+01  1.8579E+00  6.3017E+00  9.1852E+00  3.9374E+00  1.2849E+00
            -2.6525E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1646.83781513079        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  9.6955E-01  9.7973E-01  1.0668E+00  1.0147E+00  9.4549E-01  9.9321E-01  8.7752E-01  8.8851E-01  9.8962E-01  8.9428E-01
             1.0471E+00
 PARAMETER:  6.9072E-02  7.9523E-02  1.6466E-01  1.1458E-01  4.3948E-02  9.3182E-02 -3.0651E-02 -1.8207E-02  8.9566E-02 -1.1742E-02
             1.4603E-01
 GRADIENT:   3.1160E+01 -5.0282E+00  1.0852E+01 -1.7057E+01  3.3671E+00  1.3625E+00  6.2233E-01  5.5649E-01 -4.3855E+00 -1.0803E+01
            -1.1154E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1647.91236882468        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.7939E-01  9.6195E-01  1.0202E+00  1.0478E+00  9.1712E-01  9.9412E-01  7.5353E-01  6.7215E-01  1.0216E+00  9.8528E-01
             1.0482E+00
 PARAMETER:  7.9170E-02  6.1203E-02  1.1995E-01  1.4673E-01  1.3488E-02  9.4099E-02 -1.8299E-01 -2.9727E-01  1.2136E-01  8.5173E-02
             1.4710E-01
 GRADIENT:   5.2361E+01  2.0540E+01  9.3862E+00  2.7954E+01 -1.3461E+01  1.2548E+00 -5.3539E-01  2.2045E-01 -5.2373E-01  5.0052E-01
            -8.6522E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1648.74010388956        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.6706E-01  1.0395E+00  8.1615E-01  9.8248E-01  8.6192E-01  9.9110E-01  9.1064E-01  4.1673E-01  1.0175E+00  8.9385E-01
             1.0586E+00
 PARAMETER:  6.6505E-02  1.3874E-01 -1.0315E-01  8.2325E-02 -4.8597E-02  9.1059E-02  6.3925E-03 -7.7532E-01  1.1738E-01 -1.2221E-02
             1.5696E-01
 GRADIENT:   1.8963E+01  7.6379E+00 -4.7177E+00  1.4794E+01  8.1715E-01 -5.5499E-01  4.6203E-01  1.1596E+00  2.2214E+00  1.9586E+00
            -1.1915E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1648.76393232074        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.6175E-01  1.0299E+00  8.1330E-01  9.8233E-01  8.5749E-01  9.9201E-01  9.2365E-01  3.5617E-01  1.0072E+00  8.9044E-01
             1.0614E+00
 PARAMETER:  6.0995E-02  1.2942E-01 -1.0666E-01  8.2171E-02 -5.3742E-02  9.1973E-02  2.0580E-02 -9.3235E-01  1.0716E-01 -1.6041E-02
             1.5961E-01
 GRADIENT:   6.5299E+00  2.3737E+00 -2.8818E+00  5.6644E+00  8.4759E-01 -3.5771E-01  3.3582E-01  7.6510E-01  9.8502E-01  1.3055E+00
            -1.5646E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1648.77425738454        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  9.5905E-01  1.0241E+00  8.0459E-01  9.8205E-01  8.5131E-01  9.9291E-01  9.3458E-01  2.7745E-01  1.0003E+00  8.8566E-01
             1.0627E+00
 PARAMETER:  5.8187E-02  1.2378E-01 -1.1743E-01  8.1885E-02 -6.0982E-02  9.2886E-02  3.2339E-02 -1.1821E+00  1.0031E-01 -2.1423E-02
             1.6082E-01
 GRADIENT:   1.0443E-01 -3.5499E-01 -1.5655E+00  7.3339E-01  7.0045E-01 -1.5937E-01  1.9793E-01  4.2781E-01  2.7172E-01  7.7132E-01
             2.9846E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1648.82306127008        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  9.5756E-01  1.0216E+00  7.8847E-01  9.8015E-01  8.4256E-01  9.9380E-01  9.4566E-01  1.2788E-01  9.9484E-01  8.7816E-01
             1.0632E+00
 PARAMETER:  5.6628E-02  1.2141E-01 -1.3766E-01  7.9955E-02 -7.1313E-02  9.3777E-02  4.4132E-02 -1.9567E+00  9.4827E-02 -2.9923E-02
             1.6131E-01
 GRADIENT:  -3.6505E+00 -1.9521E+00 -4.2965E-01 -2.4563E+00  3.9694E-01  3.9133E-02  5.2453E-02  8.6130E-02 -2.3354E-01  2.4453E-01
             4.9231E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1649.21159409328        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      561
 NPARAMETR:  9.6861E-01  1.0340E+00  7.9671E-01  9.7496E-01  8.5305E-01  1.0011E+00  9.3482E-01  1.0000E-02  1.0025E+00  8.9022E-01
             1.0669E+00
 PARAMETER:  6.8104E-02  1.3344E-01 -1.2727E-01  7.4638E-02 -5.8938E-02  1.0111E-01  3.2600E-02 -5.0716E+00  1.0245E-01 -1.6287E-02
             1.6477E-01
 GRADIENT:  -1.4250E+01 -4.2820E+00  8.3496E-01 -6.0379E+00 -8.0155E-02 -2.9561E-01  2.4386E-03  0.0000E+00 -7.4168E-01 -1.9418E-02
             1.2758E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1649.33291985081        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      736
 NPARAMETR:  9.7456E-01  1.1439E+00  7.6755E-01  9.1002E-01  8.8754E-01  1.0016E+00  8.6162E-01  1.0000E-02  1.0663E+00  9.0074E-01
             1.0649E+00
 PARAMETER:  7.4234E-02  2.3445E-01 -1.6455E-01  5.7147E-03 -1.9303E-02  1.0165E-01 -4.8937E-02 -7.0401E+00  1.6419E-01 -4.5366E-03
             1.6292E-01
 GRADIENT:  -1.5411E+00  1.1562E+00  5.6555E-01  7.3248E-01 -6.6656E-01 -1.2412E-01 -5.0808E-02  0.0000E+00 -3.0566E-01 -2.5012E-01
             4.4923E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1649.33552807969        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      915
 NPARAMETR:  9.7528E-01  1.1602E+00  7.6249E-01  8.9897E-01  8.9324E-01  1.0020E+00  8.5117E-01  1.0000E-02  1.0783E+00  9.0402E-01
             1.0647E+00
 PARAMETER:  7.4966E-02  2.4860E-01 -1.7116E-01 -6.5071E-03 -1.2902E-02  1.0199E-01 -6.1141E-02 -7.0513E+00  1.7536E-01 -9.0676E-04
             1.6267E-01
 GRADIENT:   8.1686E-04 -3.3915E-03 -1.9444E-04 -3.9265E-03 -3.7938E-04  1.3616E-04  1.1101E-03  0.0000E+00  9.8196E-04  9.4745E-04
             4.4333E-04

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1649.33552807969        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      937
 NPARAMETR:  9.7528E-01  1.1602E+00  7.6249E-01  8.9897E-01  8.9324E-01  1.0020E+00  8.5117E-01  1.0000E-02  1.0783E+00  9.0402E-01
             1.0647E+00
 PARAMETER:  7.4966E-02  2.4860E-01 -1.7116E-01 -6.5071E-03 -1.2902E-02  1.0199E-01 -6.1141E-02 -7.0513E+00  1.7536E-01 -9.0676E-04
             1.6267E-01
 GRADIENT:   8.1686E-04 -3.3915E-03 -1.9444E-04 -3.9265E-03 -3.7938E-04  1.3616E-04  1.1101E-03  0.0000E+00  9.8196E-04  9.4745E-04
             4.4333E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      937
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0245E-05 -1.9147E-02 -3.2851E-04  7.4734E-03 -2.2675E-02
 SE:             2.9805E-02  1.8886E-02  1.5392E-04  2.5305E-02  2.3786E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9973E-01  3.1066E-01  3.2822E-02  7.6773E-01  3.4044E-01

 ETASHRINKSD(%)  1.4900E-01  3.6729E+01  9.9484E+01  1.5227E+01  2.0313E+01
 ETASHRINKVR(%)  2.9778E-01  5.9968E+01  9.9997E+01  2.8135E+01  3.6500E+01
 EBVSHRINKSD(%)  4.9150E-01  3.6321E+01  9.9520E+01  1.5434E+01  1.9414E+01
 EBVSHRINKVR(%)  9.8059E-01  5.9449E+01  9.9998E+01  2.8486E+01  3.5058E+01
 RELATIVEINF(%)  9.8816E+01  1.5050E+00  2.5634E-04  3.6260E+00  6.5045E+00
 EPSSHRINKSD(%)  4.3334E+01
 EPSSHRINKVR(%)  6.7890E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1649.3355280796875     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -914.18470151594931     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.80
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1649.336       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.75E-01  1.16E+00  7.62E-01  8.99E-01  8.93E-01  1.00E+00  8.51E-01  1.00E-02  1.08E+00  9.04E-01  1.06E+00
 


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
+       -7.69E+00  4.85E+02
 
 TH 3
+        1.18E+01  2.40E+02  4.73E+02
 
 TH 4
+       -1.19E+01  3.94E+02 -1.85E+02  8.47E+02
 
 TH 5
+       -4.00E+00 -4.12E+02 -6.03E+02  2.02E+02  1.10E+03
 
 TH 6
+        6.88E-01 -1.47E+00  2.45E+00 -2.72E+00  5.42E-01  1.91E+02
 
 TH 7
+        1.13E+00  1.22E+01  4.48E+00 -1.13E+01 -1.55E+01 -3.31E-01  4.37E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.62E+00 -2.87E+01 -2.28E+01  3.51E+01 -5.39E-01 -4.31E-01  2.63E+01  0.00E+00  9.10E+01
 
 TH10
+       -3.98E+00 -7.18E+00 -5.20E+01 -1.19E+01 -5.25E+01  3.69E-01  2.34E+01  0.00E+00  3.91E+00  9.95E+01
 
 TH11
+       -9.98E+00 -1.95E+01 -3.24E+01 -3.89E+00  6.94E+00  2.59E+00  8.37E+00  0.00E+00  8.97E+00  2.17E+01  1.89E+02
 
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
 #CPUT: Total CPU Time in Seconds,       14.409
Stop Time:
Sat Sep 18 12:12:08 CDT 2021
