/*
    Zachary Patchen
    02/29/2024
    Assignment 4
    Introduction to Data Structures and Algorithms
 */
import java.util.LinkedList;
import java.util.Queue;
import java.util.Stack;

public class Assignment {
    /**
     * The time complexity of the isBalanced method is O(n), where n is the length of the input string s.
     * This is because the method iterates through each character in the string exactly once.
     * The space complexity is also O(n) in the worst case. This is because the stack can potentially store
     * all opening brackets if they are not matched by closing brackets, leading to a space usage proportional
     * to the length of the input string. In the best case (when the input string is balanced),
     * the space complexity is O(1) as the stack remains empty.
     *
     * @param s
     * @return "YES" if input is balanced and "NO" otherwise
     */
    public static String isBalanced(String s) {
        Stack<Character> stack = new Stack<>();

        for (char bracket : s.toCharArray()) {
            if (bracket == '(' || bracket == '[' || bracket == '{') {
                stack.push(bracket);
            } else {
                if (stack.isEmpty()) {
                    return "NO";  // Closing bracket with no corresponding opening bracket
                }

                char top = stack.pop();
                if ((bracket == ')' && top != '(') ||
                        (bracket == ']' && top != '[') ||
                        (bracket == '}' && top != '{')) {
                    return "NO";  // Mismatched brackets
                }
            }
        }

        return stack.isEmpty() ? "YES" : "NO";  // Check if all opening brackets are matched
    }

    /**
     * This algorithm has a linear time complexity of O(n) and linear space complexity of O(n),
     * where n is the length of the input DNA sequence. Iterating through each nucleotide in the DNA sequence
     * takes O(n) time, where n is the length of the input DNA sequence. Each operation inside the loop
     * (pushing onto the stack) is constant time O(1). Popping characters from the stack and appending
     * them to the StringBuilder in the second loop takes O(n) time. Therefore, O(N) + O(1) + O(N) is O(N)
     * <p>
     * The space complexity is primarily determined by the space used for the stack and the StringBuilder.
     * The stack space depends on the number of 'T' nucleotides in the DNA sequence. In the worst case,
     * if all nucleotides are 'T', the space complexity is O(n). The StringBuilder space also depends
     * on the length of the input DNA sequence.
     *
     * @param dnaSequence
     * @return the transcribed RNA sequence as a String
     */
    public static String transcribeToRNA(String dnaSequence) {
        // Create a stack to store characters of the transcribed RNA sequence
        Stack<Character> rnaStack = new Stack<>();
        // Iterate through each nucleotide in the DNA sequence
        for (char nucleotide : dnaSequence.toCharArray()) {
            // If the nucleotide is 'T', push 'U' onto the stack (transcription)
            if (nucleotide == 'T') {
                rnaStack.push('U');
            } else {
                // If the nucleotide is not 'T', push the nucleotide onto the stack
                rnaStack.push(nucleotide);
            }
        }
        // Create a StringBuilder to construct the final RNA sequence
        StringBuilder rnaBuilder = new StringBuilder();
        // Pop characters from the stack and append to the StringBuilder
        while (!rnaStack.isEmpty()) {
            rnaBuilder.append(rnaStack.pop());
        }
        // Reverse the StringBuilder to get the correct order of RNA sequence
        return rnaBuilder.reverse().toString();
    }

    /**
     * Iterating through each nucleotide in the RNA sequence takes O(n) time, where n is the length of the input
     * RNA sequence. Each operation inside the loop (adding to the queue, polling from the queue, appending to
     * StringBuilder) is constant time O(1). The second loop that processes any remaining nucleotides has a
     * time complexity of O(m), where m is the number of remaining nucleotides (incomplete codons). Therefore,
     * the overall time complexity is O(n + m).
     * <p>
     * The space complexity is determined by the space used for the aminoAcidSequence StringBuilder and the codonQueue.
     * The StringBuilder space complexity is O(n) since it stores the amino acid sequence. The queue space complexity
     * is O(3) in the worst case (three nucleotides per codon), which is just O(1). Therefore,
     * the overall space complexity is O(n).
     *
     * @param rnaSequence
     * @return
     */

    public static String transcribeToAminoAcids(String rnaSequence) {
        StringBuilder aminoAcidSequence = new StringBuilder();
        Queue<Character> codonQueue = new LinkedList<>();

        for (char nucleotide : rnaSequence.toCharArray()) {
            // Ignore non-RNA nucleotides
            if (isValidRNA(nucleotide)) {
                // Add the current nucleotide to the codon queue
                codonQueue.add(nucleotide);

                // Check if the codon queue has three nucleotides
                if (codonQueue.size() == 3 ) {
                    // Extract the codon and map it to the corresponding amino acid
                    String codon = codonQueue.poll() + "" + codonQueue.poll() + codonQueue.poll();
                    char aminoAcid = mapCodonToAminoAcid(codon);

                    // Append the amino acid to the result
                    aminoAcidSequence.append(aminoAcid);
                }
            }
        }

        // Handle remaining nucleotides (incomplete codons)
        while (!codonQueue.isEmpty()) {
            aminoAcidSequence.append(".");
            codonQueue.poll();
        }

        return aminoAcidSequence.toString();
    }


    /**
     * The method iterates through each character in the input infix expression once, so the time complexity is O(n),
     * where n is the length of the input infix expression. The space complexity is determined by the stack,
     * which stores operators during the conversion process. In the worst case, the stack may store all the
     * operators from the infix expression. Therefor, the space complexity is O(n).
     * @param infix expression
     * @return postfix expression
     */
    public static String infixToPostfix(String infix) {
        // StringBuilder to store the postfix expression
        StringBuilder postfix = new StringBuilder();
        // Stack to temporarily store operators during conversion
        Stack<Character> stack = new Stack<>();
        // Iterate through each character in the infix expression
        for (char ch : infix.toCharArray()) {
            // If the character is a letter or digit, append it to the postfix expression
            if (Character.isLetterOrDigit(ch)) {
                postfix.append(ch);
            } else if (ch == '(') {
                // If the character is an opening parenthesis, push it onto the stack
                stack.push(ch);
            } else if (ch == ')') {
                // If the character is a closing parenthesis,
                // pop operators from the stack and append to the postfix expression until an opening parenthesis is encountered
                while (!stack.isEmpty() && stack.peek() != '(') {
                    postfix.append(stack.pop());
                }
                stack.pop(); // Discard the opening parenthesis
            } else if (isOperator(ch)) {
                // If the character is an operator,
                // pop operators from the stack and append to the postfix expression while their precedence is greater or equal
                while (!stack.isEmpty() && getPrecedence(ch) <= getPrecedence(stack.peek())) {
                    postfix.append(stack.pop());
                }
                // Push the current operator onto the stack
                stack.push(ch);
            }
        }
        // Pop any remaining operators from the stack and append to the postfix expression
        while (!stack.isEmpty()) {
            postfix.append(stack.pop());
        }
        // Return the final postfix expression
        return postfix.toString();
    }
    private static char mapCodonToAminoAcid(String codon) {
        switch (codon) {
            case "UUU": case "UUC": return 'F';
            case "UUA": case "UUG": case "CUU": case "CUC": case "CUA": case "CUG": return 'L';
            case "AUU": case "AUC": case "AUA": return 'I';
            case "AUG": return 'M';
            case "GUU": case "GUC": case "GUA": case "GUG": return 'V';
            case "UCU": case "UCC": case "UCA": case "UCG": return 'S';
            case "CCU": case "CCC": case "CCA": case "CCG": return 'P';
            case "ACU": case "ACC": case "ACA": case "ACG": return 'T';
            case "GCU": case "GCC": case "GCA": case "GCG": return 'A';
            case "UAU": case "UAC": return 'Y';
            case "CAA": case "CAG": return 'Q';
            case "CAU": case "CAC": return 'H';
            case "AAA": case "AAG": return 'K';
            case "GAU": case "GAC": return 'D';
            case "GAA": case "GAG": return 'E';
            case "UGU": case "UGC": return 'C';
            case "UGG": return 'W';
            case "CGU": case "CGC": case "CGA": case "CGG": return 'R';
            case "AGU": case "AGC": return 'S';
            case "AGA": case "AGG": return 'R';
            case "GGU": case "GGC": case "GGA": case "GGG": return 'G';
            case "UAA": case "UAG": case "UGA": return '.'; // Stop codons
            default: return '?'; // Unknown codon
        }
    }

    private static boolean isValidRNA(char nucleotide) {
        return nucleotide == 'A' || nucleotide == 'C' || nucleotide == 'G' || nucleotide == 'U';
    }

    // Helper function to check if the given character is an operator
    private static boolean isOperator(char c) {
        return c == '+' || c == '-' || c == '*' || c == '/' || c == '^';
    }

    // Helper function to get the precedence of the operator
    private static int getPrecedence(char operator) {
        switch (operator) {
            case '^':
                return 3;
            case '*':
            case '/':
                return 2;
            case '+':
            case '-':
                return 1;
            default:
                return -1; // For '(' and ')'
        }
    }

    public static void main(String[] args) {
        // Test isBalanced method
        System.out.println("Test isBalanced:");
        System.out.println("(({})) is balanced: " + Assignment.isBalanced("(({}))")); // YES
        System.out.println("([)] is balanced: " + Assignment.isBalanced("([)]")); // NO

        // Test transcribeToRNA method
        System.out.println("\nTest transcribeToRNA:");
        System.out.println("DNA sequence: ATCGTACGTA");
        System.out.println("Transcribed RNA sequence: " + Assignment.transcribeToRNA("ATCGTACGTA"));

        // Test transcribeToAminoAcids method
        System.out.println("\nTest transcribeToAminoAcids:");
        System.out.println("RNA sequence: AUGCCGUAA");
        System.out.println("Translated Amino Acid sequence: " + Assignment.transcribeToAminoAcids("AUGCCGUAA"));

        // Test infixToPostfix method
        System.out.println("\nTest infixToPostfix:");
        String infixExpression = "a+b*(c^d-e)^(f+g*h)-i";
        String postfixExpression = Assignment.infixToPostfix(infixExpression);
        System.out.println("Infix Expression: " + infixExpression);
        System.out.println("Postfix Expression: " + postfixExpression);
        }

}

