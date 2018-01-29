from openpathsampling.netcdfplus import StorableNamedObject

class DefaultTaskScheduler(StorableNamedObject):
    def wrap_task(self, task, *args, **kwargs):
        return task(*args, **kwargs)

    def wrap_hook(self, hook, *args, **kwargs):
        # hook wrap with fire_and_forget in dask! only care about
        # side-effects
        return hook(*args, **kwargs)
