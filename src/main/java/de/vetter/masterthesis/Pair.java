package de.vetter.masterthesis;

/**
 * Trivial pair class
 * @author David Emanuel Vetter
 *
 * @param <T1>
 * @param <T2>
 */
public class Pair<T1, T2> {
	private T1 first;
	private T2 second;
	
	public Pair(T1 first, T2 second) {
		this.first = first;
		this.second = second;
	}
	
	public void setFirst(T1 newFirst) { first = newFirst; }
	public T1 getFirst() { return first; }
	
	public void setSecond(T2 newSecond) { second = newSecond; }
	public T2 getSecond() { return second; }
}
