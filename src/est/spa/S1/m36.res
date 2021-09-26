Sat Sep 25 09:51:28 CDT 2021
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
$DATA ../../../../data/spa/S1/dat36.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m36.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1662.67510165901        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.1606E+01 -6.7701E+01 -2.9020E+01 -2.9532E+01  6.6665E+01 -3.3796E+01 -7.9703E+00  2.7754E+00  2.8612E+01 -2.0463E+01
             3.1605E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1671.16456326628        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.9930E-01  1.1254E+00  9.6759E-01  9.8041E-01  1.0079E+00  1.1059E+00  1.0849E+00  9.7875E-01  7.9428E-01  1.1232E+00
             9.8760E-01
 PARAMETER:  9.9301E-02  2.1811E-01  6.7056E-02  8.0219E-02  1.0790E-01  2.0065E-01  1.8147E-01  7.8519E-02 -1.3032E-01  2.1619E-01
             8.7520E-02
 GRADIENT:   2.3365E+01  3.0879E+01 -8.8284E+00  5.1127E+01  1.5009E+01  1.7717E+01 -6.3296E+00  3.1865E+00 -4.1507E+00  1.8906E+00
            -2.4692E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1672.47120551297        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0052E+00  1.0271E+00  8.4815E-01  1.0258E+00  8.9389E-01  1.1021E+00  1.2971E+00  6.5267E-01  6.8665E-01  9.9236E-01
             9.9284E-01
 PARAMETER:  1.0522E-01  1.2677E-01 -6.4703E-02  1.2544E-01 -1.2176E-02  1.9723E-01  3.6010E-01 -3.2669E-01 -2.7593E-01  9.2331E-02
             9.2818E-02
 GRADIENT:   3.3790E+01  2.9263E+01 -1.6770E+00  4.4238E+01  9.0093E+00  1.5215E+01  1.9243E+00  8.8838E-01 -9.3617E+00 -4.8809E+00
            -9.0086E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1673.54316682776        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9266E-01  1.0195E+00  7.5040E-01  1.0034E+00  8.4658E-01  1.0791E+00  1.2333E+00  3.3850E-01  7.4472E-01  1.0094E+00
             9.8724E-01
 PARAMETER:  9.2628E-02  1.1927E-01 -1.8714E-01  1.0341E-01 -6.6549E-02  1.7615E-01  3.0969E-01 -9.8323E-01 -1.9475E-01  1.0934E-01
             8.7154E-02
 GRADIENT:   5.5765E+00  9.2836E-01 -8.8581E+00  1.2527E+01  1.0211E+01  4.4111E+00 -8.0364E-01  9.0262E-01 -2.0188E+00  4.6040E+00
             1.4734E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1673.63791569646        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.9067E-01  1.0140E+00  7.0534E-01  9.9750E-01  8.1137E-01  1.0721E+00  1.2401E+00  2.3957E-01  7.5448E-01  9.4690E-01
             9.8463E-01
 PARAMETER:  9.0626E-02  1.1388E-01 -2.4907E-01  9.7496E-02 -1.0903E-01  1.6964E-01  3.1522E-01 -1.3289E+00 -1.8172E-01  4.5437E-02
             8.4509E-02
 GRADIENT:   4.4902E-01 -6.7675E-01 -2.9285E+00  2.7303E+00  2.6117E+00  7.5425E-01 -2.7373E-01  4.8369E-01 -1.9256E-01  1.7353E+00
             6.2514E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1673.64622700924        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  9.9049E-01  1.0133E+00  6.8425E-01  9.9418E-01  7.9668E-01  1.0709E+00  1.2426E+00  1.8013E-01  7.5588E-01  9.2147E-01
             9.8367E-01
 PARAMETER:  9.0449E-02  1.1325E-01 -2.7944E-01  9.4165E-02 -1.2730E-01  1.6848E-01  3.1723E-01 -1.6141E+00 -1.7988E-01  1.8213E-02
             8.3533E-02
 GRADIENT:  -4.2390E-01 -8.3920E-01 -1.2303E+00  2.8642E-01  7.0721E-01 -1.2198E-02 -7.9053E-02  2.8582E-01  1.4830E-01  7.5405E-01
             2.9937E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1673.66259031923        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  9.9056E-01  1.0141E+00  6.6881E-01  9.9121E-01  7.8668E-01  1.0705E+00  1.2433E+00  1.0884E-01  7.5657E-01  9.0406E-01
             9.8300E-01
 PARAMETER:  9.0515E-02  1.1401E-01 -3.0225E-01  9.1171E-02 -1.3994E-01  1.6809E-01  3.1778E-01 -2.1179E+00 -1.7897E-01 -8.6002E-04
             8.2854E-02
 GRADIENT:  -6.5653E-01 -6.8235E-01 -1.3307E-01 -9.3779E-01 -3.0174E-01 -3.3354E-01  4.1048E-02  1.0796E-01  2.4686E-01  4.9885E-02
             4.9092E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1673.77863934268        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      538
 NPARAMETR:  9.9081E-01  1.0160E+00  6.8991E-01  9.9315E-01  7.9970E-01  1.0711E+00  1.2416E+00  1.0000E-02  7.5438E-01  9.1969E-01
             9.8372E-01
 PARAMETER:  9.0771E-02  1.1584E-01 -2.7120E-01  9.3125E-02 -1.2352E-01  1.6872E-01  3.1644E-01 -4.7917E+00 -1.8186E-01  1.6282E-02
             8.3590E-02
 GRADIENT:  -4.3961E+01 -2.4399E+00  3.9572E+00 -1.1393E+01 -5.2432E+00 -1.3270E+01 -1.9233E+00  0.0000E+00 -1.8325E+00 -1.9330E+00
            -5.2194E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1674.62986116551        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      715
 NPARAMETR:  1.0118E+00  9.3070E-01  7.3696E-01  1.0532E+00  7.9381E-01  1.1004E+00  1.3495E+00  1.0000E-02  7.2971E-01  9.5925E-01
             9.8420E-01
 PARAMETER:  1.1168E-01  2.8185E-02 -2.0522E-01  1.5186E-01 -1.3091E-01  1.9570E-01  3.9976E-01 -7.3895E+00 -2.1510E-01  5.8394E-02
             8.4077E-02
 GRADIENT:  -2.2377E+00  1.3203E+00 -5.4898E-01  2.3562E+00  2.1945E-01 -8.0934E-01 -3.8910E-01  0.0000E+00 -5.0796E-01  2.6234E-01
             2.1673E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1674.64267110492        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      890
 NPARAMETR:  1.0127E+00  8.9825E-01  7.4552E-01  1.0712E+00  7.8588E-01  1.1022E+00  1.3934E+00  1.0000E-02  7.2207E-01  9.5833E-01
             9.8342E-01
 PARAMETER:  1.1265E-01 -7.3059E-03 -1.9367E-01  1.6874E-01 -1.4095E-01  1.9733E-01  4.3178E-01 -6.7827E+00 -2.2563E-01  5.7442E-02
             8.3278E-02
 GRADIENT:  -5.1480E-03  2.6874E-03  7.3204E-03 -5.7513E-03 -6.9490E-03  3.3355E-03 -1.0637E-03  0.0000E+00 -3.5726E-03 -3.0559E-03
            -1.6240E-03

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1674.64267110492        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      912
 NPARAMETR:  1.0127E+00  8.9825E-01  7.4552E-01  1.0712E+00  7.8588E-01  1.1022E+00  1.3934E+00  1.0000E-02  7.2207E-01  9.5833E-01
             9.8342E-01
 PARAMETER:  1.1265E-01 -7.3059E-03 -1.9367E-01  1.6874E-01 -1.4095E-01  1.9733E-01  4.3178E-01 -6.7827E+00 -2.2563E-01  5.7442E-02
             8.3278E-02
 GRADIENT:  -5.1480E-03  2.6874E-03  7.3204E-03 -5.7513E-03 -6.9490E-03  3.3355E-03 -1.0637E-03  0.0000E+00 -3.5726E-03 -3.0559E-03
            -1.6240E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      912
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.3460E-05  5.5874E-03 -4.6450E-04 -1.0572E-02 -7.8717E-03
 SE:             2.9869E-02  2.2371E-02  1.8726E-04  2.2698E-02  2.3966E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9964E-01  8.0278E-01  1.3121E-02  6.4136E-01  7.4257E-01

 ETASHRINKSD(%)  1.0000E-10  2.5054E+01  9.9373E+01  2.3960E+01  1.9712E+01
 ETASHRINKVR(%)  1.0000E-10  4.3830E+01  9.9996E+01  4.2179E+01  3.5539E+01
 EBVSHRINKSD(%)  3.4404E-01  2.5028E+01  9.9434E+01  2.4387E+01  1.7852E+01
 EBVSHRINKVR(%)  6.8690E-01  4.3792E+01  9.9997E+01  4.2827E+01  3.2518E+01
 RELATIVEINF(%)  9.8934E+01  4.6499E+00  3.5301E-04  4.7529E+00  6.7187E+00
 EPSSHRINKSD(%)  4.3963E+01
 EPSSHRINKVR(%)  6.8599E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1674.6426711049162     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -939.49184454117801     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.73
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.88
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1674.643       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  8.98E-01  7.46E-01  1.07E+00  7.86E-01  1.10E+00  1.39E+00  1.00E-02  7.22E-01  9.58E-01  9.83E-01
 


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
+        8.85E+02
 
 TH 2
+       -7.01E+00  4.17E+02
 
 TH 3
+        1.68E+01  1.66E+02  7.13E+02
 
 TH 4
+       -1.06E+01  4.05E+02 -3.85E+02  1.08E+03
 
 TH 5
+       -2.95E+00 -3.11E+02 -8.44E+02  4.06E+02  1.32E+03
 
 TH 6
+        5.47E-01 -2.20E+00  4.03E+00 -2.29E+00 -1.42E+00  1.62E+02
 
 TH 7
+        1.25E+00  3.31E+01 -8.28E+00 -1.52E+01  2.36E+00  3.51E-02  3.85E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.73E+00 -2.39E+01 -4.36E+01  1.60E+01  1.99E+01 -8.05E-01  2.00E+01  0.00E+00  1.41E+02
 
 TH10
+       -7.44E-01 -9.14E+00 -6.70E+01 -2.29E+01 -4.86E+01 -6.63E-01  8.36E+00  0.00E+00  1.67E+01  9.48E+01
 
 TH11
+       -4.96E+00 -1.44E+01 -3.66E+01 -1.64E+00  8.17E+00  2.19E+00  3.67E+00  0.00E+00  1.46E+01  2.13E+01  2.23E+02
 
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
 #CPUT: Total CPU Time in Seconds,       14.664
Stop Time:
Sat Sep 25 09:51:45 CDT 2021
