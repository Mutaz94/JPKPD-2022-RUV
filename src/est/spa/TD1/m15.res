Wed Sep 29 17:56:43 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat15.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m15.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1648.71588666753        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.2837E+02 -6.6249E+01 -1.5089E+01 -5.3665E+01  7.8482E+01  5.8476E+01 -1.1828E+01 -9.3626E+00 -2.1995E+01 -7.7977E+00
            -4.1655E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1662.68938424029        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      195
 NPARAMETR:  1.0481E+00  1.0538E+00  9.7297E-01  1.0591E+00  9.5921E-01  9.2450E-01  1.0499E+00  1.0640E+00  1.0734E+00  9.8625E-01
             1.1001E+00
 PARAMETER:  1.4701E-01  1.5242E-01  7.2602E-02  1.5737E-01  5.8358E-02  2.1500E-02  1.4874E-01  1.6203E-01  1.7083E-01  8.6160E-02
             1.9538E-01
 GRADIENT:   3.8533E+00 -7.9098E+00 -9.7869E+00  1.0150E+01  1.8319E+01 -3.9219E+00 -2.5238E+00 -4.1862E+00  1.5145E+00  2.3237E+00
             1.5407E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1663.73679141059        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0508E+00  1.2959E+00  9.6155E-01  9.2658E-01  1.0310E+00  9.3667E-01  9.2456E-01  1.4391E+00  1.1898E+00  9.2782E-01
             1.0932E+00
 PARAMETER:  1.4951E-01  3.5917E-01  6.0787E-02  2.3743E-02  1.3057E-01  3.4575E-02  2.1560E-02  4.6404E-01  2.7375E-01  2.5085E-02
             1.8912E-01
 GRADIENT:   8.2818E+00  2.3751E+01  1.5359E+00  2.6318E+01 -3.3535E+00  9.4867E-01 -1.9580E+00 -1.2221E+00 -3.9713E-01 -4.7393E+00
            -2.7068E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1664.30768899852        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.0464E+00  1.3143E+00  1.0536E+00  8.9546E-01  1.0855E+00  9.3410E-01  8.7019E-01  1.6131E+00  1.2293E+00  1.0176E+00
             1.0977E+00
 PARAMETER:  1.4535E-01  3.7328E-01  1.5225E-01 -1.0413E-02  1.8208E-01  3.1830E-02 -3.9045E-02  5.7816E-01  3.0648E-01  1.1742E-01
             1.9325E-01
 GRADIENT:  -7.2746E-01  1.4415E+00  1.6350E+00  1.3814E+00 -3.8674E+00  2.8709E-01 -5.2559E-01 -3.2678E-01 -3.7100E-01  2.3931E-01
             3.0300E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1664.38384793943        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      722
 NPARAMETR:  1.0476E+00  1.4711E+00  9.9479E-01  7.9712E-01  1.1324E+00  9.3426E-01  8.2091E-01  1.7456E+00  1.3314E+00  1.0452E+00
             1.0974E+00
 PARAMETER:  1.4648E-01  4.8598E-01  9.4779E-02 -1.2675E-01  2.2435E-01  3.1995E-02 -9.7342E-02  6.5708E-01  3.8622E-01  1.4420E-01
             1.9291E-01
 GRADIENT:   9.6293E-01  1.0337E+01  4.3537E+00  4.2744E+00 -7.8832E+00  1.9440E-01 -7.2272E-01 -7.1598E-01 -1.0774E+00  7.5040E-02
            -9.5258E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1664.48744535728        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      898
 NPARAMETR:  1.0486E+00  1.6238E+00  8.6737E-01  6.9557E-01  1.1644E+00  9.3424E-01  8.0264E-01  1.7954E+00  1.4471E+00  1.0565E+00
             1.0980E+00
 PARAMETER:  1.4743E-01  5.8480E-01 -4.2294E-02 -2.6302E-01  2.5223E-01  3.1978E-02 -1.1985E-01  6.8524E-01  4.6958E-01  1.5494E-01
             1.9351E-01
 GRADIENT:   2.0750E+00  1.1339E+01  3.7181E+00  4.3167E+00 -5.5257E+00 -5.8732E-02  1.5643E-01 -4.6315E-01 -7.7412E-01 -1.4898E-01
            -1.1364E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1664.58163640315        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1073
 NPARAMETR:  1.0492E+00  1.8185E+00  6.5608E-01  5.6154E-01  1.1936E+00  9.3460E-01  7.8050E-01  1.7600E+00  1.6363E+00  1.0619E+00
             1.0999E+00
 PARAMETER:  1.4806E-01  6.9802E-01 -3.2147E-01 -4.7708E-01  2.7696E-01  3.2359E-02 -1.4782E-01  6.6529E-01  5.9244E-01  1.6009E-01
             1.9519E-01
 GRADIENT:   2.2471E+00  6.7829E+00  1.5687E+00  3.6526E+00 -1.2897E+00 -2.6149E-01  7.7587E-01 -9.5253E-03 -1.7981E-01 -1.7016E-01
            -8.5276E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1664.64133897681        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1261
 NPARAMETR:  1.0492E+00  1.8105E+00  6.4851E-01  5.6028E-01  1.1907E+00  9.3522E-01  7.7594E-01  1.7312E+00  1.6373E+00  1.0615E+00
             1.1015E+00
 PARAMETER:  1.4803E-01  6.9358E-01 -3.3307E-01 -4.7933E-01  2.7451E-01  3.3023E-02 -1.5368E-01  6.4884E-01  5.9303E-01  1.5969E-01
             1.9664E-01
 GRADIENT:   2.1832E+00 -4.2676E+00  6.0347E-01  3.9000E-01  9.2923E-01  1.9883E-02 -7.4974E-02  7.9802E-02  6.9176E-02  1.7582E-01
            -6.1354E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1664.64537320353        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1442
 NPARAMETR:  1.0479E+00  1.8108E+00  6.4129E-01  5.6075E-01  1.1874E+00  9.3473E-01  7.7601E-01  1.7037E+00  1.6419E+00  1.0576E+00
             1.1004E+00
 PARAMETER:  1.4680E-01  6.9376E-01 -3.4427E-01 -4.7849E-01  2.7173E-01  3.2503E-02 -1.5359E-01  6.3280E-01  5.9586E-01  1.5600E-01
             1.9568E-01
 GRADIENT:  -1.0109E+00 -3.2564E+00  4.0075E-01  1.5311E+00  1.1845E+00 -2.2521E-01 -2.6169E-01 -5.7548E-03  7.5423E-01  2.7613E-02
            -3.5795E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1664.65703182269        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1625
 NPARAMETR:  1.0491E+00  1.8109E+00  6.3767E-01  5.5976E-01  1.1848E+00  9.3532E-01  7.7783E-01  1.6951E+00  1.6346E+00  1.0546E+00
             1.1011E+00
 PARAMETER:  1.4790E-01  6.9380E-01 -3.4994E-01 -4.8025E-01  2.6955E-01  3.3137E-02 -1.5125E-01  6.2774E-01  5.9138E-01  1.5314E-01
             1.9634E-01
 GRADIENT:   1.7460E+00 -3.5675E+00  7.4962E-01  4.3726E-01  3.1339E-01  3.0117E-02 -1.8861E-01 -2.0482E-02  5.4966E-02  3.0797E-02
            -1.0620E-01

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1664.65951571641        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:     1722
 NPARAMETR:  1.0494E+00  1.8087E+00  6.3842E-01  5.5930E-01  1.1845E+00  9.3538E-01  7.7901E-01  1.6915E+00  1.6352E+00  1.0544E+00
             1.1013E+00
 PARAMETER:  1.4824E-01  6.9260E-01 -3.4876E-01 -4.8108E-01  2.6932E-01  3.3200E-02 -1.4974E-01  6.2562E-01  5.9176E-01  1.5299E-01
             1.9651E-01
 GRADIENT:  -3.5746E-01  1.1888E+00 -1.2041E+04 -1.7168E+01  3.2740E-01  2.8570E-04  1.9971E-02  6.6600E+03 -7.1202E+03  6.3046E-02
             1.7354E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1722
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.5458E-04 -3.8936E-02 -3.2005E-02  3.1660E-02 -5.0349E-02
 SE:             2.9817E-02  2.2138E-02  1.2202E-02  2.2843E-02  2.0813E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7714E-01  7.8615E-02  8.7163E-03  1.6575E-01  1.5560E-02

 ETASHRINKSD(%)  1.0872E-01  2.5835E+01  5.9122E+01  2.3474E+01  3.0272E+01
 ETASHRINKVR(%)  2.1732E-01  4.4995E+01  8.3290E+01  4.1438E+01  5.1380E+01
 EBVSHRINKSD(%)  6.1077E-01  2.5256E+01  6.1611E+01  2.4835E+01  2.7606E+01
 EBVSHRINKVR(%)  1.2178E+00  4.4133E+01  8.5263E+01  4.3502E+01  4.7591E+01
 RELATIVEINF(%)  9.8723E+01  3.3281E+00  1.4505E+00  3.5388E+00  1.6233E+01
 EPSSHRINKSD(%)  4.5376E+01
 EPSSHRINKVR(%)  7.0162E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1664.6595157164118     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -929.50868915267358     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1664.660       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.81E+00  6.38E-01  5.59E-01  1.18E+00  9.35E-01  7.79E-01  1.69E+00  1.64E+00  1.05E+00  1.10E+00
 


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
+        1.14E+03
 
 TH 2
+       -7.57E+00  6.72E+04
 
 TH 3
+       -2.02E+02  8.98E+01  1.51E+02
 
 TH 4
+        2.51E+06  3.06E+02 -1.41E+02  7.52E+02
 
 TH 5
+       -3.77E+00 -1.29E+02  1.48E+06  1.31E+02  3.90E+02
 
 TH 6
+        5.89E-01 -1.34E+00 -2.74E+02 -2.96E+00 -9.90E-01  2.22E+02
 
 TH 7
+        6.71E-01 -5.21E+00 -8.92E+00  3.35E+06 -9.37E+00 -7.12E-01  1.14E+02
 
 TH 8
+        4.37E+01 -9.44E+00 -1.54E+01  8.48E+00  1.90E-01  5.76E+01  7.84E+00  4.51E+00
 
 TH 9
+       -4.59E+01 -1.42E+01 -4.86E+05  4.49E+01 -4.78E+00 -6.33E+01  1.27E+01  1.02E+05  1.12E+05
 
 TH10
+        3.68E-01 -7.07E+00  2.92E+06 -1.31E+02 -5.94E+01  3.22E-01  1.12E+01  3.88E+00 -3.41E+01  4.04E+06
 
 TH11
+       -8.77E+00 -1.31E+01 -2.18E+06  2.46E-01 -8.62E+00  3.52E+00  8.33E+00 -3.52E+02  3.94E+02 -3.01E+06  2.24E+06
 
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
 #CPUT: Total CPU Time in Seconds,       31.389
Stop Time:
Wed Sep 29 17:57:16 CDT 2021
