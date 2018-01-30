from openpathsampling.netcdfplus import StorableNamedObject

try:
    import dask
    import dask.distributed
except ImportError:
    HAS_DASK = False
else:
    HAS_DASK = True

class TaskScheduler(StorableNamedObject):
    def wrap_task(self, task, *args, **kwargs):
        return task(*args, **kwargs)

    def wrap_hook(self, hook, *args, **kwargs):
        # hook wrap with fire_and_forget in dask! only care about
        # side-effects
        return hook(*args, **kwargs)

class DaskTaskScheduler(TaskScheduler):
    def __init__(self, client):
        import dask.distributed  # raise ImportError if missing
        self.client = client

    def wrap_task(self, task, *args, **kwargs):
        return self.client.submit(task, *args, pure=False, **kwargs)

    def wrap_hook(self, hook, *args, **kwargs):
        return self.client.submit(hook, *args, **kwargs)

