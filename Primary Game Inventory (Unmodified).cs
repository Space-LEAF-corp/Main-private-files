public class GameInventory : IInventory
{
    private Dictionary<string, int> items = new();

    public bool AddItem(IItem item, int quantity)
    {
        if (!items.ContainsKey(item.Id))
            items[item.Id] = 0;

        // Game’s original limit logic
        if (items[item.Id] + quantity > item.MaxStack)
            return false;

        items[item.Id] += quantity;
        return true;
    }

    public bool RemoveItem(string itemId, int quantity)
    {
        if (!items.ContainsKey(itemId) || items[itemId] < quantity)
            return false;

        items[itemId] -= quantity;
        return true;
    }

    public int GetQuantity(string itemId)
    {
        return items.TryGetValue(itemId, out int qty) ? qty : 0;
    }
}
