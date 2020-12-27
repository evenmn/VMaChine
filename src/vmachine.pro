TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG += qt
LIBS += -llapack -lblas

QT += core
QT += charts
# greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

# target.path = $$/home/evenmn/VMaChine/src-eigen/

# Remove possible other optimization flags
QMAKE_CXXFLAGS_RELEASE -= -O
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2

# Add the desired -O3 if not present
QMAKE_CXXFLAGS_RELEASE *= -O3

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE       = $$QMAKE_CXX
QMAKE_CXX_DEBUG         = $$QMAKE_CXX
QMAKE_LINK              = $$QMAKE_CXX
QMAKE_CC                = mpicc

QMAKE_CFLAGS           += $$system(mpicc  --showme:compile)
QMAKE_LFLAGS           += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS         += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

# Load internal source files
SOURCES += src-eigen/main.cpp \
    src-eigen/Activation/activation.cpp \
    src-eigen/Activation/elu.cpp \
    src-eigen/Activation/leakyrelu.cpp \
    src-eigen/Activation/purelinear.cpp \
    src-eigen/Activation/relu.cpp \
    src-eigen/Activation/sigmoid.cpp \
    src-eigen/InitialWeights/customized.cpp \
    src-eigen/InitialWeights/randomnormal.cpp \
    src-eigen/InitialWeights/randomuniform.cpp \
    src-eigen/Layer/dense.cpp \
    src-eigen/Layer/input.cpp \
    src-eigen/Layer/layer.cpp \
    src-eigen/Layer/output.cpp \
    src-eigen/WaveFunctions/doubleproduct.cpp \
    src-eigen/WaveFunctions/drbmproduct.cpp \
    src-eigen/WaveFunctions/fnn.cpp \
    src-eigen/WaveFunctions/rbmproduct.cpp \
    src-eigen/system.cpp \
    src-eigen/sampler.cpp \
    src-eigen/Hamiltonians/hamiltonian.cpp \
    src-eigen/Hamiltonians/harmonicoscillator.cpp \
    src-eigen/Hamiltonians/atomicnucleus.cpp \
    src-eigen/Hamiltonians/doublewell.cpp \
    src-eigen/InitialStates/initialstate.cpp \
    src-eigen/InitialStates/randomuniform.cpp \
    src-eigen/InitialStates/randomnormal.cpp \
    src-eigen/InitialWeights/initialweights.cpp \
    src-eigen/InitialWeights/constant.cpp \
    src-eigen/InitialWeights/automatize.cpp \
    src-eigen/Metropolis/metropolis.cpp \
    src-eigen/Metropolis/bruteforce.cpp \
    src-eigen/Metropolis/importancesampling.cpp \
    src-eigen/WaveFunctions/wavefunction.cpp \
    src-eigen/WaveFunctions/gaussian.cpp \
    src-eigen/WaveFunctions/slaterdeterminant.cpp \
    src-eigen/WaveFunctions/hydrogenlike.cpp \
    src-eigen/WaveFunctions/padejastrow.cpp \
    src-eigen/WaveFunctions/simplejastrow.cpp \
    src-eigen/WaveFunctions/rbmgaussian.cpp \
    src-eigen/WaveFunctions/partlyrestricted.cpp \
    src-eigen/Optimization/optimization.cpp \
    src-eigen/Optimization/gradientdescent.cpp \
    src-eigen/Optimization/sgd.cpp \
    src-eigen/Optimization/asgd.cpp \
    src-eigen/Optimization/adam.cpp \
    src-eigen/RNG/rng.cpp \
    src-eigen/RNG/mersennetwister.cpp \
    src-eigen/Basis/basis.cpp \
    src-eigen/Basis/hermite.cpp \
    src-eigen/Basis/hydrogenorbital.cpp \
    src-eigen/Basis/none.cpp \
    src-eigen/Basis/hartreefock.cpp \
    src-eigen/Basis/hermiteexpansion.cpp

# Load external source files
SOURCES += \
    src-eigen/block/c++/blocker.cpp

# Load internal header files
HEADERS += \
    src-eigen/Activation/activation.h \
    src-eigen/Activation/elu.h \
    src-eigen/Activation/leakyrelu.h \
    src-eigen/Activation/purelinear.h \
    src-eigen/Activation/relu.h \
    src-eigen/Activation/sigmoid.h \
    src-eigen/Eigen/Cholesky \
    src-eigen/Eigen/CholmodSupport \
    src-eigen/Eigen/Core \
    src-eigen/Eigen/Dense \
    src-eigen/Eigen/Eigen \
    src-eigen/Eigen/Eigenvalues \
    src-eigen/Eigen/Geometry \
    src-eigen/Eigen/Householder \
    src-eigen/Eigen/IterativeLinearSolvers \
    src-eigen/Eigen/Jacobi \
    src-eigen/Eigen/LU \
    src-eigen/Eigen/MetisSupport \
    src-eigen/Eigen/OrderingMethods \
    src-eigen/Eigen/PaStiXSupport \
    src-eigen/Eigen/PardisoSupport \
    src-eigen/Eigen/QR \
    src-eigen/Eigen/QtAlignedMalloc \
    src-eigen/Eigen/SPQRSupport \
    src-eigen/Eigen/SVD \
    src-eigen/Eigen/Sparse \
    src-eigen/Eigen/SparseCholesky \
    src-eigen/Eigen/SparseCore \
    src-eigen/Eigen/SparseLU \
    src-eigen/Eigen/SparseQR \
    src-eigen/Eigen/StdDeque \
    src-eigen/Eigen/StdList \
    src-eigen/Eigen/StdVector \
    src-eigen/Eigen/SuperLUSupport \
    src-eigen/Eigen/UmfPackSupport \
    src-eigen/Eigen/src/Cholesky/LDLT.h \
    src-eigen/Eigen/src/Cholesky/LLT.h \
    src-eigen/Eigen/src/Cholesky/LLT_LAPACKE.h \
    src-eigen/Eigen/src/CholmodSupport/CholmodSupport.h \
    src-eigen/Eigen/src/Core/Array.h \
    src-eigen/Eigen/src/Core/ArrayBase.h \
    src-eigen/Eigen/src/Core/ArrayWrapper.h \
    src-eigen/Eigen/src/Core/Assign.h \
    src-eigen/Eigen/src/Core/AssignEvaluator.h \
    src-eigen/Eigen/src/Core/Assign_MKL.h \
    src-eigen/Eigen/src/Core/BandMatrix.h \
    src-eigen/Eigen/src/Core/Block.h \
    src-eigen/Eigen/src/Core/BooleanRedux.h \
    src-eigen/Eigen/src/Core/CommaInitializer.h \
    src-eigen/Eigen/src/Core/ConditionEstimator.h \
    src-eigen/Eigen/src/Core/CoreEvaluators.h \
    src-eigen/Eigen/src/Core/CoreIterators.h \
    src-eigen/Eigen/src/Core/CwiseBinaryOp.h \
    src-eigen/Eigen/src/Core/CwiseNullaryOp.h \
    src-eigen/Eigen/src/Core/CwiseTernaryOp.h \
    src-eigen/Eigen/src/Core/CwiseUnaryOp.h \
    src-eigen/Eigen/src/Core/CwiseUnaryView.h \
    src-eigen/Eigen/src/Core/DenseBase.h \
    src-eigen/Eigen/src/Core/DenseCoeffsBase.h \
    src-eigen/Eigen/src/Core/DenseStorage.h \
    src-eigen/Eigen/src/Core/Diagonal.h \
    src-eigen/Eigen/src/Core/DiagonalMatrix.h \
    src-eigen/Eigen/src/Core/DiagonalProduct.h \
    src-eigen/Eigen/src/Core/Dot.h \
    src-eigen/Eigen/src/Core/EigenBase.h \
    src-eigen/Eigen/src/Core/ForceAlignedAccess.h \
    src-eigen/Eigen/src/Core/Fuzzy.h \
    src-eigen/Eigen/src/Core/GeneralProduct.h \
    src-eigen/Eigen/src/Core/GenericPacketMath.h \
    src-eigen/Eigen/src/Core/GlobalFunctions.h \
    src-eigen/Eigen/src/Core/IO.h \
    src-eigen/Eigen/src/Core/Inverse.h \
    src-eigen/Eigen/src/Core/Map.h \
    src-eigen/Eigen/src/Core/MapBase.h \
    src-eigen/Eigen/src/Core/MathFunctions.h \
    src-eigen/Eigen/src/Core/MathFunctionsImpl.h \
    src-eigen/Eigen/src/Core/Matrix.h \
    src-eigen/Eigen/src/Core/MatrixBase.h \
    src-eigen/Eigen/src/Core/NestByValue.h \
    src-eigen/Eigen/src/Core/NoAlias.h \
    src-eigen/Eigen/src/Core/NumTraits.h \
    src-eigen/Eigen/src/Core/PermutationMatrix.h \
    src-eigen/Eigen/src/Core/PlainObjectBase.h \
    src-eigen/Eigen/src/Core/Product.h \
    src-eigen/Eigen/src/Core/ProductEvaluators.h \
    src-eigen/Eigen/src/Core/Random.h \
    src-eigen/Eigen/src/Core/Redux.h \
    src-eigen/Eigen/src/Core/Ref.h \
    src-eigen/Eigen/src/Core/Replicate.h \
    src-eigen/Eigen/src/Core/ReturnByValue.h \
    src-eigen/Eigen/src/Core/Reverse.h \
    src-eigen/Eigen/src/Core/Select.h \
    src-eigen/Eigen/src/Core/SelfAdjointView.h \
    src-eigen/Eigen/src/Core/SelfCwiseBinaryOp.h \
    src-eigen/Eigen/src/Core/Solve.h \
    src-eigen/Eigen/src/Core/SolveTriangular.h \
    src-eigen/Eigen/src/Core/SolverBase.h \
    src-eigen/Eigen/src/Core/StableNorm.h \
    src-eigen/Eigen/src/Core/Stride.h \
    src-eigen/Eigen/src/Core/Swap.h \
    src-eigen/Eigen/src/Core/Transpose.h \
    src-eigen/Eigen/src/Core/Transpositions.h \
    src-eigen/Eigen/src/Core/TriangularMatrix.h \
    src-eigen/Eigen/src/Core/VectorBlock.h \
    src-eigen/Eigen/src/Core/VectorwiseOp.h \
    src-eigen/Eigen/src/Core/Visitor.h \
    src-eigen/Eigen/src/Core/arch/AVX/Complex.h \
    src-eigen/Eigen/src/Core/arch/AVX/MathFunctions.h \
    src-eigen/Eigen/src/Core/arch/AVX/PacketMath.h \
    src-eigen/Eigen/src/Core/arch/AVX/TypeCasting.h \
    src-eigen/Eigen/src/Core/arch/AVX512/MathFunctions.h \
    src-eigen/Eigen/src/Core/arch/AVX512/PacketMath.h \
    src-eigen/Eigen/src/Core/arch/AltiVec/Complex.h \
    src-eigen/Eigen/src/Core/arch/AltiVec/MathFunctions.h \
    src-eigen/Eigen/src/Core/arch/AltiVec/PacketMath.h \
    src-eigen/Eigen/src/Core/arch/CUDA/Complex.h \
    src-eigen/Eigen/src/Core/arch/CUDA/Half.h \
    src-eigen/Eigen/src/Core/arch/CUDA/MathFunctions.h \
    src-eigen/Eigen/src/Core/arch/CUDA/PacketMath.h \
    src-eigen/Eigen/src/Core/arch/CUDA/PacketMathHalf.h \
    src-eigen/Eigen/src/Core/arch/CUDA/TypeCasting.h \
    src-eigen/Eigen/src/Core/arch/Default/ConjHelper.h \
    src-eigen/Eigen/src/Core/arch/Default/Settings.h \
    src-eigen/Eigen/src/Core/arch/NEON/Complex.h \
    src-eigen/Eigen/src/Core/arch/NEON/MathFunctions.h \
    src-eigen/Eigen/src/Core/arch/NEON/PacketMath.h \
    src-eigen/Eigen/src/Core/arch/SSE/Complex.h \
    src-eigen/Eigen/src/Core/arch/SSE/MathFunctions.h \
    src-eigen/Eigen/src/Core/arch/SSE/PacketMath.h \
    src-eigen/Eigen/src/Core/arch/SSE/TypeCasting.h \
    src-eigen/Eigen/src/Core/arch/ZVector/Complex.h \
    src-eigen/Eigen/src/Core/arch/ZVector/MathFunctions.h \
    src-eigen/Eigen/src/Core/arch/ZVector/PacketMath.h \
    src-eigen/Eigen/src/Core/functors/AssignmentFunctors.h \
    src-eigen/Eigen/src/Core/functors/BinaryFunctors.h \
    src-eigen/Eigen/src/Core/functors/NullaryFunctors.h \
    src-eigen/Eigen/src/Core/functors/StlFunctors.h \
    src-eigen/Eigen/src/Core/functors/TernaryFunctors.h \
    src-eigen/Eigen/src/Core/functors/UnaryFunctors.h \
    src-eigen/Eigen/src/Core/products/GeneralBlockPanelKernel.h \
    src-eigen/Eigen/src/Core/products/GeneralMatrixMatrix.h \
    src-eigen/Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h \
    src-eigen/Eigen/src/Core/products/GeneralMatrixMatrixTriangular_BLAS.h \
    src-eigen/Eigen/src/Core/products/GeneralMatrixMatrix_BLAS.h \
    src-eigen/Eigen/src/Core/products/GeneralMatrixVector.h \
    src-eigen/Eigen/src/Core/products/GeneralMatrixVector_BLAS.h \
    src-eigen/Eigen/src/Core/products/Parallelizer.h \
    src-eigen/Eigen/src/Core/products/SelfadjointMatrixMatrix.h \
    src-eigen/Eigen/src/Core/products/SelfadjointMatrixMatrix_BLAS.h \
    src-eigen/Eigen/src/Core/products/SelfadjointMatrixVector.h \
    src-eigen/Eigen/src/Core/products/SelfadjointMatrixVector_BLAS.h \
    src-eigen/Eigen/src/Core/products/SelfadjointProduct.h \
    src-eigen/Eigen/src/Core/products/SelfadjointRank2Update.h \
    src-eigen/Eigen/src/Core/products/TriangularMatrixMatrix.h \
    src-eigen/Eigen/src/Core/products/TriangularMatrixMatrix_BLAS.h \
    src-eigen/Eigen/src/Core/products/TriangularMatrixVector.h \
    src-eigen/Eigen/src/Core/products/TriangularMatrixVector_BLAS.h \
    src-eigen/Eigen/src/Core/products/TriangularSolverMatrix.h \
    src-eigen/Eigen/src/Core/products/TriangularSolverMatrix_BLAS.h \
    src-eigen/Eigen/src/Core/products/TriangularSolverVector.h \
    src-eigen/Eigen/src/Core/util/BlasUtil.h \
    src-eigen/Eigen/src/Core/util/Constants.h \
    src-eigen/Eigen/src/Core/util/DisableStupidWarnings.h \
    src-eigen/Eigen/src/Core/util/ForwardDeclarations.h \
    src-eigen/Eigen/src/Core/util/MKL_support.h \
    src-eigen/Eigen/src/Core/util/Macros.h \
    src-eigen/Eigen/src/Core/util/Memory.h \
    src-eigen/Eigen/src/Core/util/Meta.h \
    src-eigen/Eigen/src/Core/util/NonMPL2.h \
    src-eigen/Eigen/src/Core/util/ReenableStupidWarnings.h \
    src-eigen/Eigen/src/Core/util/StaticAssert.h \
    src-eigen/Eigen/src/Core/util/XprHelper.h \
    src-eigen/Eigen/src/Eigenvalues/ComplexEigenSolver.h \
    src-eigen/Eigen/src/Eigenvalues/ComplexSchur.h \
    src-eigen/Eigen/src/Eigenvalues/ComplexSchur_LAPACKE.h \
    src-eigen/Eigen/src/Eigenvalues/EigenSolver.h \
    src-eigen/Eigen/src/Eigenvalues/GeneralizedEigenSolver.h \
    src-eigen/Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h \
    src-eigen/Eigen/src/Eigenvalues/HessenbergDecomposition.h \
    src-eigen/Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h \
    src-eigen/Eigen/src/Eigenvalues/RealQZ.h \
    src-eigen/Eigen/src/Eigenvalues/RealSchur.h \
    src-eigen/Eigen/src/Eigenvalues/RealSchur_LAPACKE.h \
    src-eigen/Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h \
    src-eigen/Eigen/src/Eigenvalues/SelfAdjointEigenSolver_LAPACKE.h \
    src-eigen/Eigen/src/Eigenvalues/Tridiagonalization.h \
    src-eigen/Eigen/src/Geometry/AlignedBox.h \
    src-eigen/Eigen/src/Geometry/AngleAxis.h \
    src-eigen/Eigen/src/Geometry/EulerAngles.h \
    src-eigen/Eigen/src/Geometry/Homogeneous.h \
    src-eigen/Eigen/src/Geometry/Hyperplane.h \
    src-eigen/Eigen/src/Geometry/OrthoMethods.h \
    src-eigen/Eigen/src/Geometry/ParametrizedLine.h \
    src-eigen/Eigen/src/Geometry/Quaternion.h \
    src-eigen/Eigen/src/Geometry/Rotation2D.h \
    src-eigen/Eigen/src/Geometry/RotationBase.h \
    src-eigen/Eigen/src/Geometry/Scaling.h \
    src-eigen/Eigen/src/Geometry/Transform.h \
    src-eigen/Eigen/src/Geometry/Translation.h \
    src-eigen/Eigen/src/Geometry/Umeyama.h \
    src-eigen/Eigen/src/Geometry/arch/Geometry_SSE.h \
    src-eigen/Eigen/src/Householder/BlockHouseholder.h \
    src-eigen/Eigen/src/Householder/Householder.h \
    src-eigen/Eigen/src/Householder/HouseholderSequence.h \
    src-eigen/Eigen/src/IterativeLinearSolvers/BasicPreconditioners.h \
    src-eigen/Eigen/src/IterativeLinearSolvers/BiCGSTAB.h \
    src-eigen/Eigen/src/IterativeLinearSolvers/ConjugateGradient.h \
    src-eigen/Eigen/src/IterativeLinearSolvers/IncompleteCholesky.h \
    src-eigen/Eigen/src/IterativeLinearSolvers/IncompleteLUT.h \
    src-eigen/Eigen/src/IterativeLinearSolvers/IterativeSolverBase.h \
    src-eigen/Eigen/src/IterativeLinearSolvers/LeastSquareConjugateGradient.h \
    src-eigen/Eigen/src/IterativeLinearSolvers/SolveWithGuess.h \
    src-eigen/Eigen/src/Jacobi/Jacobi.h \
    src-eigen/Eigen/src/LU/Determinant.h \
    src-eigen/Eigen/src/LU/FullPivLU.h \
    src-eigen/Eigen/src/LU/InverseImpl.h \
    src-eigen/Eigen/src/LU/PartialPivLU.h \
    src-eigen/Eigen/src/LU/PartialPivLU_LAPACKE.h \
    src-eigen/Eigen/src/LU/arch/Inverse_SSE.h \
    src-eigen/Eigen/src/MetisSupport/MetisSupport.h \
    src-eigen/Eigen/src/OrderingMethods/Amd.h \
    src-eigen/Eigen/src/OrderingMethods/Eigen_Colamd.h \
    src-eigen/Eigen/src/OrderingMethods/Ordering.h \
    src-eigen/Eigen/src/PaStiXSupport/PaStiXSupport.h \
    src-eigen/Eigen/src/PardisoSupport/PardisoSupport.h \
    src-eigen/Eigen/src/QR/ColPivHouseholderQR.h \
    src-eigen/Eigen/src/QR/ColPivHouseholderQR_LAPACKE.h \
    src-eigen/Eigen/src/QR/CompleteOrthogonalDecomposition.h \
    src-eigen/Eigen/src/QR/FullPivHouseholderQR.h \
    src-eigen/Eigen/src/QR/HouseholderQR.h \
    src-eigen/Eigen/src/QR/HouseholderQR_LAPACKE.h \
    src-eigen/Eigen/src/SPQRSupport/SuiteSparseQRSupport.h \
    src-eigen/Eigen/src/SVD/BDCSVD.h \
    src-eigen/Eigen/src/SVD/JacobiSVD.h \
    src-eigen/Eigen/src/SVD/JacobiSVD_LAPACKE.h \
    src-eigen/Eigen/src/SVD/SVDBase.h \
    src-eigen/Eigen/src/SVD/UpperBidiagonalization.h \
    src-eigen/Eigen/src/SparseCholesky/SimplicialCholesky.h \
    src-eigen/Eigen/src/SparseCholesky/SimplicialCholesky_impl.h \
    src-eigen/Eigen/src/SparseCore/AmbiVector.h \
    src-eigen/Eigen/src/SparseCore/CompressedStorage.h \
    src-eigen/Eigen/src/SparseCore/ConservativeSparseSparseProduct.h \
    src-eigen/Eigen/src/SparseCore/MappedSparseMatrix.h \
    src-eigen/Eigen/src/SparseCore/SparseAssign.h \
    src-eigen/Eigen/src/SparseCore/SparseBlock.h \
    src-eigen/Eigen/src/SparseCore/SparseColEtree.h \
    src-eigen/Eigen/src/SparseCore/SparseCompressedBase.h \
    src-eigen/Eigen/src/SparseCore/SparseCwiseBinaryOp.h \
    src-eigen/Eigen/src/SparseCore/SparseCwiseUnaryOp.h \
    src-eigen/Eigen/src/SparseCore/SparseDenseProduct.h \
    src-eigen/Eigen/src/SparseCore/SparseDiagonalProduct.h \
    src-eigen/Eigen/src/SparseCore/SparseDot.h \
    src-eigen/Eigen/src/SparseCore/SparseFuzzy.h \
    src-eigen/Eigen/src/SparseCore/SparseMap.h \
    src-eigen/Eigen/src/SparseCore/SparseMatrix.h \
    src-eigen/Eigen/src/SparseCore/SparseMatrixBase.h \
    src-eigen/Eigen/src/SparseCore/SparsePermutation.h \
    src-eigen/Eigen/src/SparseCore/SparseProduct.h \
    src-eigen/Eigen/src/SparseCore/SparseRedux.h \
    src-eigen/Eigen/src/SparseCore/SparseRef.h \
    src-eigen/Eigen/src/SparseCore/SparseSelfAdjointView.h \
    src-eigen/Eigen/src/SparseCore/SparseSolverBase.h \
    src-eigen/Eigen/src/SparseCore/SparseSparseProductWithPruning.h \
    src-eigen/Eigen/src/SparseCore/SparseTranspose.h \
    src-eigen/Eigen/src/SparseCore/SparseTriangularView.h \
    src-eigen/Eigen/src/SparseCore/SparseUtil.h \
    src-eigen/Eigen/src/SparseCore/SparseVector.h \
    src-eigen/Eigen/src/SparseCore/SparseView.h \
    src-eigen/Eigen/src/SparseCore/TriangularSolver.h \
    src-eigen/Eigen/src/SparseLU/SparseLU.h \
    src-eigen/Eigen/src/SparseLU/SparseLUImpl.h \
    src-eigen/Eigen/src/SparseLU/SparseLU_Memory.h \
    src-eigen/Eigen/src/SparseLU/SparseLU_Structs.h \
    src-eigen/Eigen/src/SparseLU/SparseLU_SupernodalMatrix.h \
    src-eigen/Eigen/src/SparseLU/SparseLU_Utils.h \
    src-eigen/Eigen/src/SparseLU/SparseLU_column_bmod.h \
    src-eigen/Eigen/src/SparseLU/SparseLU_column_dfs.h \
    src-eigen/Eigen/src/SparseLU/SparseLU_copy_to_ucol.h \
    src-eigen/Eigen/src/SparseLU/SparseLU_gemm_kernel.h \
    src-eigen/Eigen/src/SparseLU/SparseLU_heap_relax_snode.h \
    src-eigen/Eigen/src/SparseLU/SparseLU_kernel_bmod.h \
    src-eigen/Eigen/src/SparseLU/SparseLU_panel_bmod.h \
    src-eigen/Eigen/src/SparseLU/SparseLU_panel_dfs.h \
    src-eigen/Eigen/src/SparseLU/SparseLU_pivotL.h \
    src-eigen/Eigen/src/SparseLU/SparseLU_pruneL.h \
    src-eigen/Eigen/src/SparseLU/SparseLU_relax_snode.h \
    src-eigen/Eigen/src/SparseQR/SparseQR.h \
    src-eigen/Eigen/src/StlSupport/StdDeque.h \
    src-eigen/Eigen/src/StlSupport/StdList.h \
    src-eigen/Eigen/src/StlSupport/StdVector.h \
    src-eigen/Eigen/src/StlSupport/details.h \
    src-eigen/Eigen/src/SuperLUSupport/SuperLUSupport.h \
    src-eigen/Eigen/src/UmfPackSupport/UmfPackSupport.h \
    src-eigen/Eigen/src/misc/Image.h \
    src-eigen/Eigen/src/misc/Kernel.h \
    src-eigen/Eigen/src/misc/RealSvd2x2.h \
    src-eigen/Eigen/src/misc/blas.h \
    src-eigen/Eigen/src/misc/lapack.h \
    src-eigen/Eigen/src/misc/lapacke.h \
    src-eigen/Eigen/src/misc/lapacke_mangling.h \
    src-eigen/Eigen/src/plugins/ArrayCwiseBinaryOps.h \
    src-eigen/Eigen/src/plugins/ArrayCwiseUnaryOps.h \
    src-eigen/Eigen/src/plugins/BlockMethods.h \
    src-eigen/Eigen/src/plugins/CommonCwiseBinaryOps.h \
    src-eigen/Eigen/src/plugins/CommonCwiseUnaryOps.h \
    src-eigen/Eigen/src/plugins/MatrixCwiseBinaryOps.h \
    src-eigen/Eigen/src/plugins/MatrixCwiseUnaryOps.h \
    src-eigen/InitialWeights/customized.h \
    src-eigen/InitialWeights/randomnormal.h \
    src-eigen/InitialWeights/randomuniform.h \
    src-eigen/Layer/dense.h \
    src-eigen/Layer/input.h \
    src-eigen/Layer/layer.h \
    src-eigen/Layer/output.h \
    src-eigen/WaveFunctions/doubleproduct.h \
    src-eigen/WaveFunctions/drbmproduct.h \
    src-eigen/WaveFunctions/fnn.h \
    src-eigen/WaveFunctions/rbmproduct.h \
    src-eigen/main.h \
    src-eigen/system.h \
    src-eigen/sampler.h \
    src-eigen/Hamiltonians/hamiltonian.h \
    src-eigen/Hamiltonians/harmonicoscillator.h \
    src-eigen/Hamiltonians/atomicnucleus.h \
    src-eigen/Hamiltonians/doublewell.h \
    src-eigen/InitialStates/initialstate.h \
    src-eigen/InitialStates/randomuniform.h \
    src-eigen/InitialStates/randomnormal.h \
    src-eigen/InitialWeights/initialweights.h \
    src-eigen/InitialWeights/constant.h \
    src-eigen/InitialWeights/automatize.h \
    src-eigen/Metropolis/metropolis.h \
    src-eigen/Metropolis/bruteforce.h \
    src-eigen/Metropolis/importancesampling.h \
    src-eigen/WaveFunctions/wavefunction.h \
    src-eigen/WaveFunctions/gaussian.h \
    src-eigen/WaveFunctions/slaterdeterminant.h \
    src-eigen/WaveFunctions/hydrogenlike.h \
    src-eigen/WaveFunctions/padejastrow.h \
    src-eigen/WaveFunctions/simplejastrow.h \
    src-eigen/WaveFunctions/rbmgaussian.h \
    src-eigen/WaveFunctions/partlyrestricted.h \
    src-eigen/Optimization/optimization.h \
    src-eigen/Optimization/gradientdescent.h \
    src-eigen/Optimization/sgd.h \
    src-eigen/Optimization/asgd.h \
    src-eigen/Optimization/adam.h \
    src-eigen/RNG/rng.h \
    src-eigen/RNG/mersennetwister.h \
    src-eigen/Basis/basis.h \
    src-eigen/Basis/hermite.h \
    src-eigen/Basis/hydrogenorbital.h \
    src-eigen/Basis/none.h \
    src-eigen/Basis/hartreefock.h \
    src-eigen/Basis/hermiteexpansion.h

# Load external header files
HEADERS += \
    src-eigen/block/c++/blocker.h

DISTFILES += \
    src-eigen/Eigen/CMakeLists.txt
