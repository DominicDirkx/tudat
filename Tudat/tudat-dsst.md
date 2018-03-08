# **tudat-dsst**: review of changes

This document includes a review of the changes performed in Tudat during the implementation of a propagator based on the Draper Semi-analytical Satellite Theory (DSST) [Semianalytic satellite theory, Danielson (1995)].

**This document does _NOT_ review any of the files in the directory *Tudat/Astrodynamics/Propagators/DSST/*.**

## Tudat/

### Tudat/Astrodynamics/

#### Tudat/Astrodynamics/Aerodynamics/

##### Tudat/Astrodynamics/Aerodynamics/sw19571001.txt
Updated space weather file for the NLRMSISE00 atmosphere model. Downloaded from https://celestrak.com/SpaceData/sw19571001.txt on 2 May 2017.



#### Tudat/Astrodynamics/BasicAstrodynamics/


##### Tudat/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h, .cpp

* Added functions to compute the periapsis altitude of a body given its Keplerian (Cartesian) state, and the central body's average radius (and gravitational parameter).
> For central bodies that cannot be modeled as a perfect sphere, the results may be inaccurate. In this case, it would be necessary to add an implementation in which the sub-satellite point of the propagated body at periapsis passage and the shape model of the central body are used to determine the radius at that point precisely.


##### Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h, .cpp

* Added functions to convert between modified equinoctial elements (currently implemented in `tudat-master`) and equinoctial elements (new in `tudat-dsst`).


##### Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h

* Added functions to convert between Cartesian and Keplerian elements and equinoctial elements, as defined in [Semianalytic satellite theory, Danielson (1995), §2.1.5].

* Added functions to convert between true, eccentric and mean longitudes.


##### Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h

* Added indeces for the new equinoctial elements, as defined in [Semianalytic satellite theory, Danielson (1995), §2.1.5].

* Added an enumeration for the type of the fast variable (mean, eccentric or true).
> In semi-analytical satellite theory, some perturbations can be modeled more easily by using a certain type of fast variable (either mean, eccentric or true). Thus, functions taking states provided in equinoctial elements as input, also have an optional input to specify the type of the fast variable used in the provided state. By default, the type of the fast variable is mean, as this the one specified in the definition of the equinoctial element set in [Semianalytic satellite theory, Danielson (1995), §2.1.5].



#### Tudat/Astrodynamics/ElectroMagnetism/

##### Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureAcceleration.h

* The `RadiationPressureInterface` used to create a `CannonBallRadiationPressureAcceleration` can be now stored as a member of the latter.
> In the constructors of `CannonBallRadiationPressureAcceleration`, `radiationPressureInterface` is the last argument and has a default value of `NULL`, ensuring backwards compatibility.

* Added a function to get the current distance between the accelerated body and the source.

* Added getters for some of the members of `CannonBallRadiationPressureAcceleration`.
> These functions have not been documented.



#### Tudat/Astrodynamics/Gravitation/

##### Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModelBase.h

* Added a getter for the value of the current gravitational parameter.




#### Tudat/Astrodynamics/Propagators/

##### Tudat/Astrodynamics/Propagators/CMakeLists.txt

* Added header and source files in the directory *Tudat/Astrodynamics/Propagators/DSST/* to the `tudat_propagators` library.

* Added DSST unit tests.


##### Tudat/Astrodynamics/Propagators/dynamicsStateDerivativeModel.h

* Added a boolean indicating whether the termination condition should also be evaluated at each of the intermediate substeps performed by the associated integrator or only after each integration step has been completed.
> By default, the value of this member is `false`, ensuring backwards compatibility.


##### Tudat/Astrodynamics/Propagators/integrateEquations.h

* Added a `propagationTerminationReason` so that it is possible to output whether the propagation was terminated because the termination condition was reached or because and error was thrown (and caught) during propagation.
> In the method `integrateEquations`, the `propagationTerminationReason` is passed (and returned) by reference, so this is code-breaking. However, all the calls to `integrateEquations` in Tudat have been updated. Existing Tudat applications in which the function `integrateEquations` is called explicitly won't compile.

* Added the boolean `assessPropagationTerminationConditionDuringIntegrationSubsteps` as an input to the method `integrateEquations`.
> This is the last argument of the method and has a default value of `false`, ensuring backwards compatibility.

* Updated the method `integrateEquationsFromIntegrator` to support:
  * Updating of the `propagationTerminationReason`.
  * Stopping the propagation when the integrator's `getPropagationShouldTerminate()` method returns `true`, which can only happen if the integrator has been configured to check the propagation termination condition after each sub-step.


##### Tudat/Astrodynamics/Propagators/nBodyDSSTStateDerivative.h

* New file defining the class for the DSST propagator.


##### Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h

* Added `dsst` to the `TranslationalPropagatorType` enumeration.



#### Tudat/Basics/

##### Tudat/Basics/utilities.h

* Added a function returning the result of subtracting the values of two `boost::function`'s.



### Tudat/InputOutput/

##### Tudat/InputOutput/CMakeLists.txt

* Added `mapTextFileReader` to `tudat_input_output` library and corresponding unit tests.


##### Tudat/InputOutput/mapTextFileReader.h, .cpp

* New files defining a function that returns a map by reading the values from a text file. The user can specify the type of the keys `KeyType` and the scalar type of the values `ScalarValueType`. The returned object will be of type `std::map< KeyType, std::vector< ScalarValueType > >`.
> This function can read files such as:
```
1	1	2
2	1	2	4
3	1	2	4	8
4	1	2	4	8	16
```
where the first value of each row is used as the key and the rest of elements are stored as a vector of the specified `ScalarValueType`.


##### Tudat/InputOutput/UnitTests/unitTestMapTextFileReader.cpp

* Unit test for the functionality defined in `Tudat/InputOutput/mapTextFileReader.cpp`.


### Tudat/Mathematics/

#### Tudat/Mathematics/NumericalIntegrators/

##### Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h

* Updated to add support for stopping the propagation if the termination condition is reached while evaluating one of the intermediate sub-steps necessary to integrate to the next step.
> A `propagationTerminationFunction_` has been added, which will always be evaluated after each sub-step, and the result of this evaluation will be stored in `propagationShouldTerminate_`. By default, the termination function is `boost::lambda::constant( false )`, but it is possible to use the method `setPropagationTerminationFunction` to set it equal to the termination conditions provided by the user. When the condition is reached, further sub-steps are not computed and the propagator can use the method `getPropagationShouldTerminate( )` to discard the last result and stop the propagation.


##### Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h

* Updated to add support for evaluating the propagation termination condition after each intermediate sub-step (k1, k2, k3, k4).


##### Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h

* Support for evaluating the propagation termination condition after each intermediate sub-step has **NOT** been implemented yet.
> In this integrator, steps can be rejected if the error is not within bounds. The termination condition may be reached during one of the evaluations of the state derivative within an integration step, but then that step may be rejected if the error is too large. If the computation of further sub-steps was interrupted, it would not be possible to determine wether the step is within acceptable error. Not interrupting the process may lead to undefined state for some propagators (DSST).


#### Tudat/Mathematics/NumericalQuadrature/

##### Tudat/Mathematics/NumericalQuadrature/CMakeLists.txt

* Added the Gaussian quadrature to the `tudat_numerical_quadrature` library and corresponding unit tests.


##### Tudat/Mathematics/NumericalQuadrature/UnitTests/unitTestGaussianQuadrature.cpp

* Added unit test for the Gaussian quadrature.


##### Tudat/Mathematics/NumericalQuadrature/gaussianQuadrature.h

* Class for Gaussian quadrature.
> It is not possible for the user to provide custom input files for the Gaussian nodes and weights, so the quadrature is currently limited to order 64 (64 nodes).


##### Tudat/Mathematics/NumericalQuadrature/gaussianNodes.txt

* List of *unique* nodes for the Gaussian quadrature up to order 64.

> When the order is an odd number, the node 0 is always considered. This is not included in this file, but generated by the `GaussianQuadrature` class.

> Additionally, for all the nodes in this file, there is another node –the conjugate– with the same absolute value but different sign. This is not included in this file, but generated by the `GaussianQuadrature` class.


##### Tudat/Mathematics/NumericalQuadrature/gaussianWeights.txt

* List of *unique* weights for the Gaussian quadrature up to order 64.

> For **even** order: the nth weight factor is used for the nth node provided in the file `gaussianNodes.txt` and its conjugate.

> For **odd** order: the first weight factor is used for node 0; the nth weight factor (for n > 1) is used for the (n-1)th node provided in the file `gaussianNodes.txt` and its conjugate.



### Tudat/SimulationSetup/

#### Tudat/SimulationSetup/EnvironmentSetup/

##### Tudat/SimulationSetup/EnvironmentSetup/createAtmosphereModel.cpp

* Now, if the atmosphere settings have type `nrlmsise00`, they can be of the base class `AtmosphereSettings` –the default space weather file "sw19571001.txt" will be used– or of the new class `NRLMSISE00AtmosphereSettings`, through which a custom space weather file provided by the user can be specified.


##### Tudat/SimulationSetup/EnvironmentSetup/createAtmosphereModel.h

* Added the class `NRLMSISE00AtmosphereSettings`, which has a `spaceWeatherFile_` member.


#### Tudat/SimulationSetup/PropagationSetup/

##### Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.cpp

* In the method `createCannonballRadiationPressureAcceleratioModel`, the `radiationPressureInterface` is passed as input.

##### Tudat/SimulationSetup/PropagationSetup/createStateDerivativeModel.h

* The `case dsst` has been in `createTranslationalStateDerivativeModel`.


##### Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h

* If a DSST propagator is being used, `dynamicsStateDerivative_->assessPropagationTerminationConditionDuringIntegrationSubsteps` is set to `true`.

* Support for the `PropagationTerminationReason`.


##### Tudat/SimulationSetup/PropagationSetup/propagationOutput.h, .cpp

* Added implementations for the new `PropagationDependentVariables`:
  * `periapsis_altitude_dependent_variable`.
  * `dsst_mean_element_rates`.
  * `dsst_short_period_terms`.
  * `dsst_computation_times`.

* Added function `evaluateTrivariateFunction`, which evaluates three `boost::function`'s and uses the values as input for another function.

* Previously, scalar dependent variables were outputted first, followed by all the vectorial dependent variables. Now all the dependent variables are outputted in the order in which they are added to the `dependentVariables_` container of the `DependentVariableSaveSettings`.
> To achieve this, scalar variables are transformed to a vectorial variable of size one, and thus all the variables can be treated in the same way.


##### Tudat/SimulationSetup/PropagationSetup/propagationOutputSettings.h, .cpp

* Added four new dependent variables and their names.

* Added class `ForceIdentifiableDependentVariableSaveSettings`, used to request the `dsst_mean_element_rates`, `dsst_short_period_terms` and/or `dsst_computation_times` only for a specific force (e.g. Earth drag, SRP, etc.).
> `ForceIdentifier` is defined (currently in the namespace of the DSST propagator) as a struct containing the name of the body causing the acceleration (`std::string`) and the type of the acceleration (`basic_astrodynamics::AvailableAcceleration`).


##### Tudat/SimulationSetup/PropagationSetup/propagationSettings.h

* Added class `DSSTTranslationalStatePropagatorSettings` to provide custom settings for the DSST propagator, such as number of elements in the series expansions for third-body attraction and SRP, number of nodes for the Gaussian quadrature, etc.


##### Tudat/SimulationSetup/PropagationSetup/propagationTermination.h

* Added `PropagationTerminationReason`:
  * `unknown_reason`.
  * `termination_condition_reached`.
  * `runtime_error_caught`.


##### Tudat/SimulationSetup/PropagationSetup/variationalEquationsSolver.h

* Added support for the `PropagationTerminationReason`.
