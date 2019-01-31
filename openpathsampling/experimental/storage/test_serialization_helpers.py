from collections import namedtuple
import json
from .serialization_helpers import *
from .serialization_helpers import _uuids_from_table_row
import numpy as np
import pytest

def toy_uuid_maker(name):
    return int(hash(name))

def uuid_encode(name):
    return "UUID(" + str(toy_uuid_maker(name)) + ")"

class MockUUIDObject(object):
    attr_list = ['name', 'normal_attr', 'obj_attr', 'list_attr',
                 'dict_attr', 'lazy_attr']
    schema = [('dict_attr', 'uuid'), ('list_attr', 'list_uuid'),
              ('obj_attr', 'uuid'), ('lazy_attr', 'lazy'),
              ('normal_attr', 'str')]
    def __init__(self, name, normal_attr=None, obj_attr=None,
                 list_attr=None, dict_attr=None, lazy_attr=None):
        self.name = name
        self.__uuid__ = toy_uuid_maker(name)
        self.dict_attr = dict_attr
        self.list_attr = list_attr
        self.obj_attr = obj_attr
        self.normal_attr = normal_attr
        self.lazy_attr = lazy_attr

    def to_dict(self):
        return {
            'name': self.name,
            'obj_attr': self.obj_attr,
            'list_attr': self.list_attr,
            'dict_attr': self.dict_attr,
            'normal_attr': self.normal_attr,
            'lazy_attr': self.lazy_attr
        }

    @classmethod
    def from_dict(cls, dct):
        # set UUID after
        return cls(name=None, **dct)

def create_test_objects():
    obj_int = MockUUIDObject(name='int', normal_attr=5)
    obj_str = MockUUIDObject(name='str', normal_attr='foo')
    obj_np = MockUUIDObject(name='np', normal_attr=np.array([1.0, 2.0]))
    obj_obj = MockUUIDObject(name='obj', obj_attr=obj_int)
    obj_lst = MockUUIDObject(name='lst', list_attr=[obj_int, obj_str])
    obj_dct = MockUUIDObject(name='dct', dict_attr={'foo': obj_str,
                                                    obj_int: obj_np})
    obj_nest = MockUUIDObject(
        name='nest',
        dict_attr={'bar': [obj_str, {obj_int: [obj_np, obj_obj]}]}
    )
    obj_repeat = MockUUIDObject('rep', list_attr=[obj_int, [obj_int]])
    all_objects = {
        obj.name : obj
        for obj in [obj_int, obj_str, obj_np, obj_obj, obj_lst, obj_dct,
                    obj_nest, obj_repeat]
    }
    return all_objects

all_objects = create_test_objects()

class MockBackend(object):
    def _table_data_for_object(self, obj, table_name, **kwargs):
        schema_entries = self.schema[table_name]
        row_type = self.row_types[table_name]
        uuid = get_uuid(obj)
        table_idx = self.table_names.index(table_name)
        idx = len(self.tables[table_idx])
        row_kwargs = {'uuid': uuid, 'idx': idx}
        row_kwargs.update(kwargs)
        row = row_type(**row_kwargs)
        uuid_row = self.row_types['uuids'](uuid=uuid, table=table_idx,
                                           idx=idx)
        return uuid_row, row

    def __init__(self):
        self.schema = {
            'objs': [('obj_attr', 'uuid')],
            'ints': [('normal_attr', 'int')],
            'sims': [('json', 'json_obj'), ('class_idx', 'int')]
        }
        self.row_types = {
            'uuids':  namedtuple('UUIDsRow', ['uuid', 'table', 'idx']),
            'objs': namedtuple('ObjRow', ['uuid', 'idx', 'obj_attr']),
            'ints': namedtuple('IntRow', ['uuid', 'idx', 'normal_attr']),
            'sims': namedtuple('SimRow',
                               ['uuid', 'idx', 'json', 'class_idx']),
            'clss': namedtuple('ClsRow', ['idx', 'module', 'cls'])
        }
        self.obj_to_table = {
            'int': 'ints',
            'obj': 'objs',
            'dct': 'sims',
            'str': 'sims'
        }
        self.table_names = ['uuids', 'clss', 'sims', 'objs', 'ints']

        self.uuid_table = {}
        self.tables = [[] for table in self.table_names]
        uuid_row, table_row = self._table_data_for_object(
            obj=all_objects['int'],
            table_name='ints',
            normal_attr=5
        )
        self.uuid_table[uuid_row.uuid] = uuid_row
        self.tables[uuid_row.table].append(table_row)

    def load_uuids_table(self, new_uuids):
        return [self.uuid_table[uuid] for uuid in new_uuids]

    def load_table_data(self, uuid_rows):
        return [self.tables[row.table][row.idx] for row in uuid_rows]

    def uuid_row_to_table_name(self, row):
        return self.table_names[row.table]


class TestMockBackend(object):
    def setup(self):
        self.backend = MockBackend()

    @pytest.mark.parametrize(('obj_name', 'table_idx', 'idx'),
                             [('int', 4, 0)])
    def test_load_uuids_table(self, obj_name, table_idx, idx):
        # TODO: add more examples to test
        uuid = get_uuid(all_objects[obj_name])
        assert self.backend.load_uuids_table([uuid]) == \
                [(uuid, table_idx, idx)]

    def test_load_table_data(self):
        pass

    def test_uuid_row_to_table_name(self):
        pass


@pytest.mark.parametrize('obj', list(all_objects.values()))
def test_has_uuid(obj):
    assert has_uuid(obj)

def test_has_uuid_no_uuid():
    assert not has_uuid(10)
    assert not has_uuid('foo')

@pytest.mark.parametrize('name,obj', list(all_objects.items()))
def test_get_uuid(name, obj):
    assert get_uuid(obj) == str(int(hash(name)))

def test_get_uuid_none():
    assert get_uuid(None) is None

def test_get_uuid_error():
    with pytest.raises(AttributeError):
        get_uuid(10)

def test_set_uuid():
    obj = MockUUIDObject(name="test", normal_attr=10)
    fake_uuid = 100
    assert get_uuid(obj) != str(fake_uuid)
    set_uuid(obj, fake_uuid)
    assert get_uuid(obj) == str(fake_uuid)

def test_encode_uuid():
    obj = MockUUIDObject(name="test", normal_attr=10)
    set_uuid(obj, 100)
    assert encode_uuid("UUID(100)")

@pytest.mark.parametrize('obj,included_objs', [
    (all_objects['int'], []),
    (all_objects['str'], []),
    (all_objects['np'], []),
    (all_objects['obj'], [all_objects['int']]),
    (all_objects['lst'], [all_objects['int'], all_objects['str']]),
    (all_objects['dct'], [all_objects['str'], all_objects['int'],
                          all_objects['np']]),
    (all_objects['nest'], [all_objects['str'], all_objects['int'],
                           all_objects['np'], all_objects['obj']]),
    (all_objects['rep'], [all_objects['int']])
])
def test_get_all_uuids(obj, included_objs):
    expected = {str(o.__uuid__): o for o in included_objs}
    expected.update({str(obj.__uuid__): obj})
    assert get_all_uuids(obj) == expected

@pytest.mark.parametrize('obj,included_objs', [
    (all_objects['int'], []),
    (all_objects['str'], []),
    (all_objects['np'], []),
    (all_objects['obj'], []),
    (all_objects['lst'], [all_objects['str']]),
    (all_objects['dct'], [all_objects['str'], all_objects['np']]),
    (all_objects['nest'], [all_objects['str'], all_objects['np'],
                           all_objects['obj']]),
    (all_objects['rep'], [])
])
def test_get_all_uuids_with_known(obj, included_objs):
    expected = {str(o.__uuid__): o for o in included_objs}
    known_obj = all_objects['int']
    known_uuids = {str(known_obj.__uuid__): known_obj}
    if obj is not all_objects['int']:
        expected.update({str(obj.__uuid__): obj})
    assert get_all_uuids(obj, known_uuids=known_uuids) == expected

@pytest.mark.parametrize('obj,replace_dct', [
    (all_objects['int'], {'name': 'int', 'normal_attr': 5}),
    (all_objects['str'], {'name': 'str', 'normal_attr': 'foo'}),
    (all_objects['obj'], {'name': 'obj',
                          'obj_attr': uuid_encode('int')}),
    (all_objects['lst'], {'name': 'lst',
                          'list_attr': [uuid_encode('int'),
                                        uuid_encode('str')]}),
    (all_objects['dct'], {'name': 'dct',
                          'dict_attr': {
                              'foo': uuid_encode('str'),
                              uuid_encode('int'): uuid_encode('np')
                          }}),
    (all_objects['nest'], {
        'name': 'nest',
        'dict_attr': {
            'bar': [uuid_encode('str'),
                    {uuid_encode('int'): [uuid_encode('np'),
                                          uuid_encode('obj')]}]
        }}),
    (all_objects['rep'], {'name': 'rep',
                          'list_attr': [uuid_encode('int'),
                                        [uuid_encode('int')]]})
])
def test_replace_uuid(obj, replace_dct):
    after_replacement = {key: None for key in MockUUIDObject.attr_list}
    after_replacement.update(replace_dct)
    encoding = lambda x: "UUID(" + str(x) + ")"
    assert replace_uuid(obj.to_dict(), encoding) == after_replacement

def test_replace_uuid_ndarray():
    # can't use assert == with arrays
    after_replacement = {'name': 'np', 'dict_attr': None, 'lazy_attr': None,
                         'list_attr': None, 'obj_attr': None,
                         'normal_attr': np.array([1.0, 2.0])}
    encoding = lambda x: "UUID(" + str(x) + ")"
    result = replace_uuid(all_objects['np'].to_dict(), encoding)
    assert set(result.keys()) == set(after_replacement.keys())
    for key in after_replacement:
        if key == 'normal_attr':
            assert np.allclose(result[key], after_replacement[key])
        else:
            assert result[key] == after_replacement[key]

@pytest.fixture
def cache_list():
    make_cache = lambda keys: {get_uuid(all_objects[key]): all_objects[key]
                               for key in keys}
    cache_1 = make_cache(['int', 'str'])
    cache_2 = make_cache(['obj'])
    return [cache_1, cache_2]

def test_search_caches(cache_list):
    for key in ['int', 'str', 'obj']:
        uuid = get_uuid(all_objects[key])
        assert search_caches(uuid, cache_list) == all_objects[key]

    # seearch in a single cache (listify)
    for key in ['int', 'str']:
        uuid = get_uuid(all_objects[key])
        assert search_caches(uuid, cache_list[0]) == all_objects[key]

def test_search_caches_missing(cache_list):
    assert search_caches("foo", cache_list, raise_error=False) is None
    uuid = get_uuid(all_objects['obj'])
    assert search_caches(uuid, cache_list[0], raise_error=False) is None

def test_seach_caches_missing_error(cache_list):
    with pytest.raises(KeyError):
        search_caches("foo", cache_list)

    with pytest.raises(KeyError):
        uuid = get_uuid(all_objects['obj'])
        search_caches(uuid, cache_list[0])

def test_from_dict_with_uuids(cache_list):
    # this one only uses lst
    # this matches test_replace_uuid
    lst_dict = {key: None for key in MockUUIDObject.attr_list}
    lst_dict.update({'name': 'lst', 'list_attr': [uuid_encode('int'),
                                                  uuid_encode('str')]})

    lst_obj = from_dict_with_uuids(lst_dict, cache_list)
    expected = all_objects['lst']
    assert lst_obj == expected.to_dict()

def test_uuids_from_table_row():
    TableRow = namedtuple("TableRow",
                          ['uuid', 'idx'] + MockUUIDObject.attr_list)
    row = TableRow(name="row",
                   dict_attr=None,
                   list_attr=json.dumps([uuid_encode('int'),
                                         uuid_encode('str')]),
                   obj_attr=str(toy_uuid_maker('obj')),
                   lazy_attr=str(toy_uuid_maker('nest')),
                   normal_attr="bar",
                   uuid=str(toy_uuid_maker('row')),
                   idx=1)
    entries = MockUUIDObject.schema
    uuids, lazy, deps = _uuids_from_table_row(row, entries)

    assert lazy == {str(toy_uuid_maker('nest'))}
    # TODO: do we want to allow None in the UUID list? Comes from dict_attr
    assert set(uuids) == {str(toy_uuid_maker('int')),
                          str(toy_uuid_maker('str')),
                          str(toy_uuid_maker('obj')), None}
    assert len(deps) == 1
    assert set(deps[str(toy_uuid_maker('row'))]) == \
            {str(toy_uuid_maker('int')), str(toy_uuid_maker('str')),
             str(toy_uuid_maker('obj')), str(toy_uuid_maker('nest'))}

def test_schema_find_uuids():
    test_objs = create_test_objects()
    test_objs['lazy'] = test_objs['obj']
    schemas = {
        'int': [('normal_attr', 'int')],
        'str': [('normal_attr', 'str')],
        'lst': [('list_attr', 'list_uuid')],
        'obj': [('obj_attr', 'uuid')],
        'lazy': [('obj_attr', 'lazy')]
    }
    expected_newobjs = {
        'int': [],
        'str': [],
        'lst': [test_objs['int'], test_objs['str']],
        'obj': [test_objs['int']],
        'lazy': [test_objs['int']]
    }
    for (key, schema_entries) in schemas.items():
        obj = test_objs[key]
        schema_find_uuids = SchemaFindUUIDs(schema_entries)
        uuids, new_objs = schema_find_uuids(test_objs[key], cache_list=[])
        assert uuids == {get_uuid(obj): obj}
        assert new_objs == expected_newobjs[key]

def test_get_all_uuids_loading():
    pytest.skip()

def test_get_all_uuids_loading_with_existing():
    pytest.skip()

def test_dependency_dag():
    # check that we have the expected nodes, edges
    pytest.skip()

def test_dependency_dag_with_initial_dag():
    pytest.skip()

def test_get_reload_order():
    # check order, including no-dep orders
    pytest.skip()
