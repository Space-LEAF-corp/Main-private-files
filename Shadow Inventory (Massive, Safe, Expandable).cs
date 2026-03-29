public class ShadowInventory : IInventory
{
    private Dictionary<string, long> storage = new();

    public bool AddItem(IItem item, int quantity)
    {
        if (!storage.ContainsKey(item.Id))
            storage[item.Id] = 0;

        // No stack limit — this is the point
        storage[item.Id] += quantity;
        return true;
    }

    public bool RemoveItem(string itemId, int quantity)
    {
        if (!storage.ContainsKey(itemId) || storage[itemId] < quantity)
            return false;

        storage[itemId] -= quantity;
        return true;
    }

    public int GetQuantity(string itemId)
    {
        return storage.TryGetValue(itemId, out long qty) ? (int)qty : 0;
    }
}
