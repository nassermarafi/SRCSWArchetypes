from __future__ import absolute_import
__author__ = 'marafi'


def SolutionAlgorithim(OData, Dt, Tol, Steps):
    #Insert within the While loop, make sure parameter "ok" is defined
    import OpenSeesAPI
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Lower Dt: %f and Tol: %f ... "'%(Dt,Tol)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Newton Line Search ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch(Tolerance=0.8))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze %d %f ]'%(Steps,Dt)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Newton with Initial Tangent ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.Newton(Initial=True))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze %d %f ]'%(Steps,Dt)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Broyden ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.Broyden(8))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze %d %f ]'%(Steps,Dt)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying KrylovNewton ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.KrylovNewton())
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze %d %f ]'%(Steps,Dt)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

def SolutionAlgorithimV2(OData, Dt, Tol, Steps):
    #Insert within the While loop, make sure parameter "ok" is defined
    import OpenSeesAPI
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Lower Dt: %f and Tol: %f ... "'%(Dt,Tol)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Krylov... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.KrylovNewton(MaxDim = 6))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze %d %f ]'%(Steps,Dt)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying NewtonLineSearch... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch(Tolerance=0.8))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze %d %f ]'%(Steps,Dt)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying NewtonLineSearch Bisection... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch('Bisection'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze %d %f ]'%(Steps,Dt)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying NewtonLineSearch Secant... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch('Secant'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze %d %f ]'%(Steps,Dt)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying NewtonLineSearch RegulaFalsi... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch('RegulaFalsi'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze %d %f ]'%(Steps,Dt)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

def SolutionAlgorithimKrylovOnly(OData, Dt, Tol, Steps, MaxDim = 6):
    #Insert within the While loop, make sure parameter "ok" is defined
    import OpenSeesAPI
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Lower Dt: %e and Tol: %e ... "'%(Dt,Tol)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Krylov... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tol, 1000, 2))
    # OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.KrylovNewton(MaxDim = MaxDim))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze %d %e ]'%(Steps,Dt)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

def SenSolutionAlgorithim(OData, Dt, Steps, Tol = 1e-12, KrylovMaxDim = 12, MinDt = 1e-12, NoOfIterations=3000):
    import OpenSeesAPI
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set conv_tol %e'%Tol))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set max_iter %d;'%NoOfIterations))
    OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tol, 3000, 0))
    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('test EnergyIncr $conv_tol $max_iter;'))
    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('algorithm Newton;'))
    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('integrator Newmark 0.5 0.25;'))
    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('analysis Transient;'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set dt %e;'%Dt))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set min_dt %e;'%MinDt))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set n_steps %d;'%Steps))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set cur_step 1;'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set div 10.0;'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set tol 1.0e-12;'))
    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('set eigenvalue [eigen 9];'))
    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('modalDamping 0.02;'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('while {$cur_step < $n_steps} {'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tol, NoOfIterations, 0))
    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('	test EnergyIncr $conv_tol $max_iter;'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('	algorithm Newton;'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('	set ok [analyze 1 $dt];'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('	if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('		set dt_temp [expr $dt];'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('		puts "> analysis failed to converge at step $cur_step";'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('		puts "> trying KrylovNewton";'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('		algorithm KrylovNewton -maxDim %d;'%KrylovMaxDim))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('		set ok [analyze 1 $dt];'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('		if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('			set t 0.0;'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('			set mini_t 0.0;'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('			set dt_temp [expr round($dt/$div/$tol)*$tol];'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('			set mini_dt_temp 0.0;'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('			while {$t < $dt} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('				if {$dt_temp < $min_dt} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('					puts "<< model did not converge (reason: time step less than $min_dt)";'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('					puts "<< exiting safely";'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('					wipe;'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('					exit;'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('				};'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('				if {$dt_temp < [expr $dt/pow($div, 2)]} {'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tol*10, NoOfIterations, 0))
    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('					test EnergyIncr [expr $conv_tol*10.0] $max_iter;'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('				};'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('				set ok [analyze 1 $dt_temp];'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('				if {$ok == 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('					set t [expr round(($t + $dt_temp)/$tol)*$tol];'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('					set mini_t [expr round(($mini_t + $dt_temp)/$tol)*$tol];'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('					if {$mini_t >= $mini_dt_temp} {set dt_temp [expr round($dt_temp*$div/$tol)*$tol]};'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('				} else {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('					set mini_t 0.0;'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('					set mini_dt_temp [expr round($dt_temp/$tol)*$tol];'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('					set dt_temp [expr round($dt_temp/$div/$tol)*$tol];'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('				};'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('			};'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('		};'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('	};'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('	if {$cur_step % 1 == 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('		puts "Running Tim History Step: $cur_step out of %d (Sen Algo.)";'%Steps))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('	};'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('	incr cur_step;'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('};'))


def PushOverSolutionAlgorithim(OData, StepSize, Tol, ControlNode):
    #Insert within the While loop, make sure parameter "ok" is defined
    import OpenSeesAPI
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Smaller Step: %f and Tol: %f ... "'%(StepSize,Tol)))

    OData.AddObject(OpenSeesAPI.Analysis.Integrator.Static.DisplacementControl(ControlNode, 1, StepSize))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying KrylovNewton ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.KrylovNewton())
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Newton Line Search ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch(Tolerance=0.8))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Newton with Initial Tangent ... "'))
    # OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tol,1000,0))
    # OData.AddObject(OpenSeesAPI.Analysis.Algorithm.Newton(Initial=True))
    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))
    #
    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Broyden ... "'))
    # OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,1000,0))
    # OData.AddObject(OpenSeesAPI.Analysis.Algorithm.Broyden(8))
    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    # OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Newton Line Search BiSection ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch('Bisection'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Newton Line Search Secant... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch('Secant'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Newton Line Search RegulaFalsi ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch('RegulaFalsi'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

def PushOverSolutionAlgorithimDispIncr(OData, StepSize, Tol, ControlNode):
    #Insert within the While loop, make sure parameter "ok" is defined
    import OpenSeesAPI
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Smaller Step: %f and Tol: %f ... "'%(StepSize,Tol)))

    OData.AddObject(OpenSeesAPI.Analysis.Integrator.Static.DisplacementControl(ControlNode, 1, StepSize))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying KrylovNewton ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.KrylovNewton())
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Newton Line Search ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch(Tolerance=0.8))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Newton Line Search BiSection ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch('Bisection'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Newton Line Search Secant... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch('Secant'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Newton Line Search RegulaFalsi ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch('RegulaFalsi'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

def PushOverSolutionAlgorithimConstantAlgorithm(OData, StepSize, Tol, ControlNode, Iter=1000):
    import OpenSeesAPI
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Smaller Step: %f and Tol: %f ... "'%(StepSize,Tol)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.Analysis.Integrator.Static.DisplacementControl(ControlNode, 1, StepSize))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying KrylovNewton ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,1000,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.KrylovNewton())
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

def PushOverSolutionAlgorithimConstantAlgorithmDispIncr(OData, StepSize, Tol, ControlNode, NoOfIterations=1000):
    import OpenSeesAPI
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Smaller Step: %f and Tol: %f ... "'%(StepSize,Tol)))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.Analysis.Integrator.Static.DisplacementControl(ControlNode, 1, StepSize))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying KrylovNewton ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.NormDispIncr(Tol,NoOfIterations,2))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.KrylovNewton())
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

def PushOverSolutionAlgorithimConstantTol(OData, Tol, Iter=1000):
    import OpenSeesAPI
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying KrylovNewton ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,Iter,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.KrylovNewton())
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Newton Line Search ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,Iter,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch(Tolerance=0.8))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Newton Line Search BiSection ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,Iter,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch('Bisection'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Newton Line Search Secant... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,Iter,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch('Secant'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))

    OData.AddObject(OpenSeesAPI.TCL.TCLScript('if {$ok != 0} {'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('puts "Trying Newton Line Search RegulaFalsi ... "'))
    OData.AddObject(OpenSeesAPI.Analysis.Test.EnergyIncr(Tol,Iter,0))
    OData.AddObject(OpenSeesAPI.Analysis.Algorithm.NewtonLineSearch('RegulaFalsi'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('set ok [analyze 1]'))
    OData.AddObject(OpenSeesAPI.TCL.TCLScript('}'))