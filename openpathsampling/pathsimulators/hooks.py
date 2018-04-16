import openpathsampling as paths
#from .path_simulator import PathSimulator
#from openpathsampling.netcdfplus import StorableNamedObject

"""
Hooks to change class:`.PathSimulator` behavior.

These hooks group several methods together for use as part of a
:class:`.PathSimulator` ``run`` method. They allow for additional
calculations or output at several points in the simulation.
"""


class PathSimulatorHook(object):
    """Superclass for PathSimulator hooks.

    This implementation is a do-nothing hook. Subclasses should subclass the
    relevant method in order to add hooks PathSimulator objects.
    """
    implemented_for = ['before_simulation', 'before_step', 'after_step',
                       'after_simulation']

    def before_simulation(self, sim):
        pass

    def before_step(self, sim, step_number, step_info, state):
        pass

    def after_step(self, sim, step_number, step_info, state, results,
                   hook_state):
        # TODO: shouldnt this at least return hook_state,
        # such that the next hook gets a hook_state and not None?
        pass

    def after_simulation(self, sim):
        pass


class StorageHook(PathSimulatorHook):
    """Standard hook for storage.
    """
    implemented_for = ['before_simulation', 'after_step',
                       'after_simulation']

    def __init__(self, storage=None, frequency=None):
        self.storage = storage
        self.frequency = frequency

    def before_simulation(self, sim):
        if self.storage is None:
            self.storage = sim.storage
        if self.frequency is None:
            self.frequency = sim.save_frequency

    def after_step(self, sim, step_number, step_info, state, results,
                   hook_state):
        if self.storage is not None:
            self.storage.save(results)
            if step_number % self.frequency == 0:
                self.storage.sync_all()

    def after_simulation(self, sim):
        if self.storage is not None:
            sim.storage.sync_all()


class ShootFromSnapshotsOutputHook(PathSimulatorHook):
    """Default (serial) output for ShootFromSnapshotsSimulation objects.

    Updates every time a new snapshot is shot from.

    Parameters
    ----------
    output_stream : stream
        where to write the results; default ``None`` uses the simulation's
        ``output_stream``
    allow_refresh : bool
        whether to allow refresh (see :meth:`.refresh_output`); default
        ``None`` uses the simulation's value
    """
    implemented_for = ['before_simulation', 'before_step']

    def __init__(self, output_stream=None, allow_refresh=None):
        self.output_stream = output_stream
        self.allow_refresh = allow_refresh

    def before_simulation(self, sim):
        if self.output_stream is None:
            self.output_stream = sim.output_stream
        if self.allow_refresh is None:
            self.allow_refresh = sim.allow_refresh

    def before_step(self, sim, step_number, step_info, state):
        snap_num, n_snapshots, step, n_per_snapshot = step_info
        paths.tools.refresh_output(
            "Working on snapshot %d / %d; shot %d / %d" % (
                snap_num+1, n_snapshots, step+1, n_per_snapshot
            ),
            output_stream=self.output_stream,
            refresh=self.allow_refresh
        )


class LiveVisualizerHook(PathSimulatorHook):
    """
    LiveVisualizerhook for simulation objects.

    Parameters
    ----------
    live_visualizer : openpathsampling.StepVisualizer2D; default `None` uses
        the simulations live_visualizer if any
    status_update_frequency : int
        number of steps between two refreshs of the visualization;
        default `None` uses the simulations value (default=1 for PathSampling)
    """
    implemented_for = ['before_simulation', 'after_step']

    def __init__(self, live_visualizer=None, status_update_frequency=None):
        self.live_visualizer = live_visualizer
        self.status_update_frequency = status_update_frequency

    def before_simulation(self, sim):
        if self.live_visualizer is None:
            self.live_visualizer = sim.live_visualizer
        if self.status_update_frequency is None:
            self.status_update_frequency = sim.status_update_frequency

    def after_step(self, sim, step_number, step_info, state, results,
                   hook_state):
        if step_number % self.status_update_frequency == 0:
            # do we visualize this step?
            if self.live_visualizer is not None and results is not None:
                # do we visualize at all?
                self.live_visualizer.draw_ipynb(results)


class PathSamplingOutputHook(PathSimulatorHook):
    """
    Default (serial) output for PathSamplingSimulation objects.

    Updates every time a new MCStep is undertaken.

    Parameters
    ----------
    output_stream : stream
        where to write the results; default ``None`` uses the simulation's
        ``output_stream``
    allow_refresh : bool
        whether to allow refresh (see :meth:`.refresh_output`); default
        ``None`` uses the simulation's value
    """
    implemented_for = ['before_simulation', 'before_step', 'after_simulation']

    def __init__(self, output_stream=None, allow_refresh=None):
        self.output_stream = output_stream
        self.allow_refresh = allow_refresh

    def before_simulation(self, sim):
        if self.output_stream is None:
            self.output_stream = sim.output_stream
        if self.allow_refresh is None:
            self.allow_refresh = sim.allow_refresh

    def before_step(self, sim, step_number, step_info, state):
        nn, n_steps, elapsed = step_info
        paths.tools.refresh_output(
            "Working on Monte Carlo cycle number " + str(step_number)
            + "\n" + paths.tools.progress_string(nn, n_steps, elapsed),
            refresh=self.allow_refresh,
            output_stream=self.output_stream
        )

    def after_simulation(self, sim):
        paths.tools.refresh_output(
            "DONE! Completed " + str(sim.step) + " Monte Carlo cycles.\n",
            refresh=False,
            output_stream=self.output_stream
        )
