struct UnsafeView{T, N} <: DenseArray{T, N}
    dim::NTuple{N, Int}
    ptr::Ptr{T}
end

const UnsafeVectorView{T} = UnsafeView{T,1}
const UnsafeMatrixView{T} = UnsafeView{T,2}

@inline Base.size(v::UnsafeView) = v.dim
@inline Base.size(v::UnsafeVectorView, idx::Int) = idx == 1 ? v.dim[idx] : 1
@inline Base.size(v::UnsafeMatrixView, idx::Int) = idx <= 2 ? v.dim[idx] : 1
@inline Base.getindex(v::UnsafeView, idx) = unsafe_load(v.ptr, idx)
@inline Base.setindex!(v::UnsafeView, value, idx) = unsafe_store!(v.ptr, value, idx)
@inline Base.length(v::UnsafeView) = prod(v.dim)
@inline Base.pointer(v::UnsafeView) = v.ptr
@inline Base.unsafe_convert(::Type{Ptr{T}}, v::UnsafeView{T}) where {T} = v.ptr
@compat Base.IndexStyle(::Type{V}) where {V <: UnsafeView} = IndexLinear()

# For now only the kind of views needed in the package, maybe switch to ArrayViews.jl if that becomes faster.
@inline unsafe_view(parent::DenseVector{T}, range::UnitRange) where {T} = UnsafeView{T,1}((length(range),), pointer(parent) + sizeof(eltype(parent)) * (start(range) - 1))
@inline unsafe_view(parent::DenseMatrix{T}, range::UnitRange, column::Int) where {T} = UnsafeView{T,1}((length(range),), pointer(parent) + sizeof(eltype(parent)) * ((column - 1) * size(parent, 1) + start(range) - 1))
@inline unsafe_view(parent::DenseMatrix{T}, ::Colon, column::Int) where {T} = UnsafeView{T,1}((size(parent, 1),), pointer(parent) + sizeof(eltype(parent)) * (column - 1) * size(parent, 1))
@inline unsafe_view(parent::DenseMatrix{T}, ::Colon, range::UnitRange) where {T} = UnsafeView{T,2}((size(parent, 1), length(range)), pointer(parent) + sizeof(eltype(parent)) * (first(range) - 1) * size(parent, 1))
