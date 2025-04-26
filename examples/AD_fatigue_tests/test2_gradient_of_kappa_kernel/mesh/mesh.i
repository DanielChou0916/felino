[Mesh]
  [gen]
    type = FileMeshGenerator
    file = mesh.msh
  []
  [fix_node]
      type = BoundingBoxNodeSetGenerator
      input = gen
      bottom_left = '-0.1 -0.1 0'
      top_right = '0.0001 0.0001 0'
      new_boundary = 'fix_node'
  []
  #[./scale]
  #  type = TransformGenerator
  #  input = fix_node
  #  transform = SCALE
  #  vector_value ='1e-3 1e-3 1e-3'
  #[]
[]