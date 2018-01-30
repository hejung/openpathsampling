from __future__ import absolute_import
from nose.tools import (assert_equal, assert_not_equal, assert_almost_equal,
                        raises, assert_true)
from nose.plugins.skip import Skip, SkipTest

from .test_helpers import (
    true_func, assert_equal_array_array, make_1d_traj, data_filename,
    assert_items_equal
)

import openpathsampling as paths

from openpathsampling.task_schedulers import *

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

def mock_task(argument, kwargument):
    return argument + kwargument

class MockHook(object):
    def __init__(self):
        self.val = None

    def hook(self, argument, kwargument):
        self.val = argument + kwargument

class TestTaskScheduler(object):
    def setup(self):
        self.sched = TaskScheduler()

    def test_wrap_task(self):
        assert_equal(self.sched.wrap_task(mock_task, 1, kwargument=2), 3)

    def test_wrap_hook(self):
        assert_equal(self.sched.wrap_hook(mock_task, 1, kwargument=2), 3)

class TestDaskTaskScheduler(object):
    def setup(self):
        if not HAS_DASK:  # imported in the *-import of task_scheduler
            raise SkipTest
        client = dask.distributed.Client()
        self.sched = DaskTaskScheduler(client)

    def test_wrap_task(self):
        wrapped = self.sched.wrap_task(mock_task, 1, kwargument=2)
        assert_true(isinstance(wrapped, dask.distributed.Future))
        assert_equal(wrapped.result(), 3)

    def test_wrap_hook(self):
        # TODO: modify hooks to allow stored state as output?
        hook = MockHook()
        self.sched.wrap_hook(hook.hook, 1, kwargument=2)

