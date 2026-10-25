"""Reuse #261's accepted input declaration without running its caloric check."""
import ast
from pathlib import Path

def declaration(root):
    source=root/'engine/tests/equilibrium/test_trace_reliability.py'
    tree=ast.parse(source.read_text())
    items=[]
    for node in tree.body:
        if isinstance(node,ast.Import):
            if all(x.name in ('json','math') for x in node.names): items.append(node)
        elif isinstance(node,ast.ImportFrom):
            if node.module in ('pathlib','typing','epcsaft','epcsaft.equilibrium'): items.append(node)
        elif isinstance(node,ast.Assign) and isinstance(node.targets[0],ast.Name) and node.targets[0].id=='_MEA':
            items.append(node)
        elif isinstance(node,ast.FunctionDef) and node.name=='_mea_mapping':
            items.append(node)
    function=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='test_nominally_unloaded_mea_resolves_trace_ion_stationarity')
    body=[]
    for node in function.body:
        if isinstance(node,ast.FunctionDef) and node.name=='actions': break
        if isinstance(node,ast.FunctionDef) and node.name=='solve':
            node.body=[x for x in node.body if not isinstance(x,(ast.Return,ast.Assert))
                       and not (isinstance(x,ast.Assign) and isinstance(x.targets[0],ast.Name) and x.targets[0].id=='result')]
            node.body.append(ast.parse('return problem').body[0])
        body.append(node)
    body.append(ast.parse('return model, solve(353.15)').body[0])
    function.name='declaration';function.args=ast.arguments(posonlyargs=[],args=[],kwonlyargs=[],kw_defaults=[],defaults=[])
    function.body=body;function.decorator_list=[];function.returns=None
    items.append(function)
    namespace={'__file__':str(source)}
    exec(compile(ast.fix_missing_locations(ast.Module(body=items,type_ignores=[])),str(source),'exec'),namespace)
    return namespace['declaration']()
