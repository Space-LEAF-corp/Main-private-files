public interface IItem
{
    string Id { get; }
    string Name { get; }
    int MaxStack { get; }
}

public interface IInventory
{
    bool AddItem(IItem item, int quantity);
    bool RemoveItem(string itemId, int quantity);
    int GetQuantity(string itemId);
}
