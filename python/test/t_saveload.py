#! /usr/bin/env python

import openturns as ot
import otmeshing
import os
import inspect

ot.TESTPREAMBLE()
ot.Log.Show(ot.Log.NONE)

# find all instantiable classes
persistentClasses = {}
for name, obj in inspect.getmembers(otmeshing):
    if inspect.isclass(obj) and issubclass(obj, ot.PersistentObject):
        persistentClasses[obj.__name__] = obj

# save / load
fileName = "studyStd.xml"
failed = []
for cname, class_ in persistentClasses.items():
    study = ot.Study()
    study.setStorageManager(ot.XMLStorageManager(fileName))
    print(cname)
    try:
        instance = class_()
        study.add(cname, instance)
        study.save()
        study = ot.Study()
        study.setStorageManager(ot.XMLStorageManager(fileName))
        study.load()
        os.remove(fileName)
        instance = class_()
        study.fillObject(cname, instance)
        print(cname, "OK")
    except Exception as exc:
        failed += [cname]
        print("--", cname, exc)
print(f"==== {len(failed)} failures / {len(persistentClasses)} classes ====")
print(f"failed={failed}")
assert len(failed) == 0, f"{len(failed)} serialization failures: {failed}"
