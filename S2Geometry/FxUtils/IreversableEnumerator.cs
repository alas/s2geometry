namespace S2Geometry;

public interface IReversableEnumerator<T> : IEnumerator<T>
{
    bool MovePrevious();
    void SetPosition(int position);
    bool Done();
}
