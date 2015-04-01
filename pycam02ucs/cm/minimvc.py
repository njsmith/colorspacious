# This file is part of pycam02ucs
# Copyright (C) 2014 Nathaniel Smith <njs@pobox.com>
# See file LICENSE.txt for license information.

import threading
from weakref import WeakSet

# This should be used from only one thread
class Subscription(object):
    def __init__(self):
        self._event = threading.Event()

    def _fire(self):
        self._event.set()

    def wait(self, timeout=None):
        self._event.wait(timeout)
        self._event.clear()

class Trigger(object):
    def __init__(self):
        self._subscriptions = WeakSet()
        self._callbacks = set()

    def add_callback(self, f):
        self._callbacks.add(f)
        # Always call it immediately -- this is always legal (b/c as soon as
        # you call add_callback you have to be prepared for changes to
        # happen), saves having to explicitly call refresh methods in every
        # __init__, and flushes out bugs.
        f()

    def remove_callback(self, f):
        self._callbacks.remove(f)

    def add_subscription(self):
        s = Subscription()
        self._subscriptions.add(s)
        return s

    def fire(self):
        for s in self._subscriptions:
            s._fire()
        for f in self._callbacks:
            f()
